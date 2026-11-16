####################################
############## ColabDesign functions
####################################
### Import dependencies
import gc, os, re, shutil, math, pickle
import time
import warnings
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import jax
import jax.numpy as jnp
from scipy.special import softmax
from colabdesign import mk_afdesign_model, clear_mem
from colabdesign.mpnn import mk_mpnn_model
from colabdesign.af.alphafold.common import residue_constants
from colabdesign.af.loss import get_ptm, mask_loss, get_dgram_bins, _get_con_loss
from colabdesign.shared.utils import copy_dict
from .biopython_utils import hotspot_residues, calculate_clash_score, calc_ss_percentage, calculate_percentages
from .pyrosetta_utils import pr_relax, align_pdbs
from .generic_utils import update_failures
from .logging_utils import vprint
from .ipsae_utils import calculate_ipsae

def save_design_animation(af_model, output_path):
    """Render optional ColabDesign output without losing a completed design."""
    try:
        # Pinned ColabDesign uses the legacy API removed in Matplotlib 3.9.
        # pyplot.get_cmap retains its name/Colormap/None and LUT semantics.
        if not hasattr(cm, "get_cmap"):
            cm.get_cmap = plt.get_cmap
        animation = af_model.animate(dpi=150)
        with open(output_path, 'w') as stream:
            stream.write(animation)
        return True
    except Exception as error:
        warnings.warn(f"Could not save design animation {output_path}: {error}. "
                      "Continuing with relaxation and scoring.", RuntimeWarning)
        return False
    finally:
        plt.close('all')

class TrajectoryDesignRuntime:
    """AF2 weights and JIT functions for one active exact binder length."""

    def __init__(self):
        self.model = None
        self.signature = None
        self.lengths = set()
        self.initial_opt = None

    def invalidate(self):
        """Retire a runner only after its outputs and validators are finished.

        Do not call ColabDesign's clear_mem here: it deletes *all* live backend
        buffers, including buffers held by unrelated callers. Dropping the
        model and clearing JAX caches lets unreachable buffers and executables
        go without invalidating objects that still have owners. Native host and
        device allocators may retain the freed memory for reuse.
        """
        self.model = self.initial_opt = self.signature = None
        self.lengths.clear()
        # ColabDesign's bound methods and callbacks form reference cycles.
        gc.collect()
        jax.clear_caches()
        gc.collect()

    def release_outputs(self):
        """Called after the driver has consumed the finished trajectory."""
        self.model._tmp.clear()
        del self.model.aux
        self.model._inputs.pop("prev", None)

    def prepare(self, starting_pdb, chain, hotspot, length, seed, helicity, settings):
        backend = jax.default_backend()
        starting_pdb = os.fspath(starting_pdb)
        target = starting_pdb  # ColabDesign also accepts PDB text.
        if "\n" not in starting_pdb:
            stat = os.stat(starting_pdb)
            target = (os.path.realpath(starting_pdb), stat.st_size, stat.st_mtime_ns)
        signature = (backend, target, chain, hotspot, settings["af_params_dir"], settings["use_multimer_design"],
                     settings["num_recycles_design"], settings["use_rg_loss"],
                     settings["use_i_ptm_loss"], settings["use_termini_distance_loss"])
        if self.model is not None:
            # Stock CUDA prep chooses its >384-residue subbatch regime only on construction.
            crosses_boundary = (backend != "metal" and
                                (sum(self.model._lengths) > 384) != (self.model._target_len + length > 384))
            if signature != self.signature or length not in self.lengths or crosses_boundary:
                self.invalidate()

        fresh = self.model is None
        if length in self.lengths:
            vprint(f"[AF2] Reusing design runtime for binder length {length}")
        try:
            if fresh:
                self.model = mk_afdesign_model(protocol="binder", debug=False, data_dir=settings["af_params_dir"],
                    use_multimer=settings["use_multimer_design"], num_recycles=settings["num_recycles_design"],
                    best_metric="loss")
                self.initial_opt = copy_dict(self.model.opt)
            model = self.model
            # prep_inputs snapshots opt as restart defaults and runs Metal's memory preflight.
            model.opt = copy_dict(self.initial_opt)
            model.prep_inputs(pdb_filename=starting_pdb, chain=chain, binder_len=length, hotspot=hotspot,
                              seed=seed, rm_aa=settings["omit_AAs"], rm_target_seq=settings["rm_template_seq_design"],
                              rm_target_sc=settings["rm_template_sc_design"])
            model.opt["weights"].update(pae=settings["weights_pae_intra"], plddt=settings["weights_plddt"],
                                        i_pae=settings["weights_pae_inter"], con=settings["weights_con_intra"],
                                        i_con=settings["weights_con_inter"])
            model.opt["con"].update(num=settings["intra_contact_number"], cutoff=settings["intra_contact_distance"],
                                    binary=False, seqsep=9)
            model.opt["i_con"].update(num=settings["inter_contact_number"], cutoff=settings["inter_contact_distance"],
                                      binary=False)

            losses = []
            if settings["use_rg_loss"]:
                losses.append(("rg", add_rg_loss, settings["weights_rg"]))
            if settings["use_i_ptm_loss"]:
                losses.append(("i_ptm", add_i_ptm_loss, settings["weights_iptm"]))
            if settings["use_termini_distance_loss"]:
                losses.append(("NC", add_termini_distance_loss, settings["weights_termini_loss"]))
            losses.append(("helix", add_helix_loss, helicity))
            for name, register, weight in losses:
                if fresh:
                    register(model, weight)
                else:
                    model.opt["weights"][name] = weight
        except Exception:
            # Preparation mutates the model incrementally. Retire a partially prepared runner.
            self.invalidate()
            raise

        self.signature = signature
        self.lengths.add(length)
        return model


# hallucinate a binder
def binder_hallucination(design_name, starting_pdb, chain, target_hotspot_residues, length, seed, helicity_value, design_models, advanced_settings, design_paths, failure_csv, runtime=None):
    model_pdb_path = os.path.join(design_paths["Trajectory"], design_name+".pdb")
    t0 = time.time()
    if runtime is None:
        clear_mem()
        runtime = TrajectoryDesignRuntime()
    af_model = runtime.prepare(starting_pdb, chain, target_hotspot_residues or None, length, seed,
                               helicity_value, advanced_settings)
    vprint(f"[AF2] Prepared binder model and inputs in {time.time()-t0:.2f}s")

    # calculate the number of mutations to do based on the length of the protein
    greedy_tries = math.ceil(length * (advanced_settings["greedy_percentage"] / 100))

    ### start design algorithm based on selection
    if advanced_settings["design_algorithm"] == '2stage':
        # uses gradient descend to get a PSSM profile and then uses PSSM to bias the sampling of random mutations to decrease loss
        af_model.design_pssm_semigreedy(soft_iters=advanced_settings["soft_iterations"], hard_iters=advanced_settings["greedy_iterations"], tries=greedy_tries, models=design_models, 
                                        num_models=1, sample_models=advanced_settings["sample_models"], ramp_models=False, save_best=True)

    elif advanced_settings["design_algorithm"] == '3stage':
        # 3 stage design using logits, softmax, and one hot encoding
        af_model.design_3stage(soft_iters=advanced_settings["soft_iterations"], temp_iters=advanced_settings["temporary_iterations"], hard_iters=advanced_settings["hard_iterations"], 
                                num_models=1, models=design_models, sample_models=advanced_settings["sample_models"], save_best=True)

    elif advanced_settings["design_algorithm"] == 'greedy':
        # design by using random mutations that decrease loss
        af_model.design_semigreedy(advanced_settings["greedy_iterations"], tries=greedy_tries, num_models=1, models=design_models,
                                sample_models=advanced_settings["sample_models"], save_best=True)

    elif advanced_settings["design_algorithm"] == 'mcmc':
        # design by using random mutations that decrease loss
        half_life = round(advanced_settings["greedy_iterations"] / 5, 0)
        t_mcmc = 0.01
        af_model._design_mcmc(advanced_settings["greedy_iterations"], half_life=half_life, T_init=t_mcmc, mutation_rate=greedy_tries, num_models=1, models=design_models,
                                sample_models=advanced_settings["sample_models"], save_best=True)

    elif advanced_settings["design_algorithm"] == '4stage':
        # initial logits to prescreen trajectory
        print("Stage 1: Test Logits")
        t_logits = time.time()
        af_model.design_logits(iters=50, e_soft=0.9, models=design_models, num_models=1, sample_models=advanced_settings["sample_models"], save_best=True)
        vprint(f"[AF2] Stage 1 logits in {time.time()-t_logits:.2f}s")

        # determine pLDDT of best iteration according to lowest 'loss' value
        initial_plddt = get_best_plddt(af_model, length)
        
        # if best iteration has high enough confidence then continue
        if initial_plddt > 0.65:
            print("Initial trajectory pLDDT good, continuing: "+str(initial_plddt))
            soft_iterations = advanced_settings["soft_iterations"]
            temporary_iterations = advanced_settings["temporary_iterations"]
            if advanced_settings["optimise_beta"]:
                # temporarily dump model to assess secondary structure
                af_model.save_pdb(model_pdb_path)
                t_dssp = time.time()
                _, beta, *_ = calc_ss_percentage(model_pdb_path, advanced_settings, 'B')
                os.remove(model_pdb_path)
                vprint(f"[AF2] Beta assessment via DSSP in {time.time()-t_dssp:.2f}s")

                # if beta sheeted trajectory is detected then choose to optimise
                if float(beta) > 15:
                    soft_iterations += advanced_settings["optimise_beta_extra_soft"]
                    temporary_iterations += advanced_settings["optimise_beta_extra_temp"]
                    af_model.set_opt(num_recycles=advanced_settings["optimise_beta_recycles_design"])
                    print("Beta sheeted trajectory detected, optimising settings")

            # how many logit iterations left
            logits_iter = soft_iterations - 50
            if logits_iter > 0:
                print("Stage 1: Additional Logits Optimisation")
                af_model.clear_best()
                t_logits2 = time.time()
                af_model.design_logits(iters=logits_iter, e_soft=1, models=design_models, num_models=1, sample_models=advanced_settings["sample_models"],
                                    ramp_recycles=False, save_best=True)
                vprint(f"[AF2] Additional logits in {time.time()-t_logits2:.2f}s")
                af_model._tmp["seq_logits"] = af_model.aux["seq"]["logits"]
                logit_plddt = get_best_plddt(af_model, length)
                print("Optimised logit trajectory pLDDT: "+str(logit_plddt))
            else:
                logit_plddt = initial_plddt

            # perform softmax trajectory design
            if temporary_iterations > 0:
                print("Stage 2: Softmax Optimisation")
                af_model.clear_best()
                t_soft = time.time()
                af_model.design_soft(temporary_iterations, e_temp=1e-2, models=design_models, num_models=1,
                                    sample_models=advanced_settings["sample_models"], ramp_recycles=False, save_best=True)
                vprint(f"[AF2] Softmax stage in {time.time()-t_soft:.2f}s")
                softmax_plddt = get_best_plddt(af_model, length)
            else:
                softmax_plddt = logit_plddt

            # perform one hot encoding
            if softmax_plddt > 0.65:
                print("Softmax trajectory pLDDT good, continuing: "+str(softmax_plddt))
                if advanced_settings["hard_iterations"] > 0:
                    af_model.clear_best()
                    print("Stage 3: One-hot Optimisation")
                    t_hard = time.time()
                    af_model.design_hard(advanced_settings["hard_iterations"], temp=1e-2, models=design_models, num_models=1,
                                    sample_models=advanced_settings["sample_models"], dropout=False, ramp_recycles=False, save_best=True)
                    vprint(f"[AF2] One-hot stage in {time.time()-t_hard:.2f}s")
                    onehot_plddt = get_best_plddt(af_model, length)

                if onehot_plddt > 0.65:
                    # perform greedy mutation optimisation
                    print("One-hot trajectory pLDDT good, continuing: "+str(onehot_plddt))
                    if advanced_settings["greedy_iterations"] > 0:
                        print("Stage 4: PSSM Semigreedy Optimisation")
                        t_pssm = time.time()
                        af_model.design_pssm_semigreedy(soft_iters=0, hard_iters=advanced_settings["greedy_iterations"], tries=greedy_tries, models=design_models, 
                                                        num_models=1, sample_models=advanced_settings["sample_models"], ramp_models=False, save_best=True)
                        vprint(f"[AF2] PSSM semigreedy in {time.time()-t_pssm:.2f}s")

                else:
                    update_failures(failure_csv, 'Trajectory_one-hot_pLDDT')
                    print("One-hot trajectory pLDDT too low to continue: "+str(onehot_plddt))

            else:
                update_failures(failure_csv, 'Trajectory_softmax_pLDDT')
                print("Softmax trajectory pLDDT too low to continue: "+str(softmax_plddt))

        else:
            update_failures(failure_csv, 'Trajectory_logits_pLDDT')
            print("Initial trajectory pLDDT too low to continue: "+str(initial_plddt))

    else:
        print("ERROR: No valid design model selected")
        exit()
        return

    ### save trajectory PDB
    final_plddt = get_best_plddt(af_model, length)
    t_save = time.time()
    af_model.save_pdb(model_pdb_path)
    vprint(f"[AF2] Saved trajectory PDB in {time.time()-t_save:.2f}s")
    af_model.aux["log"]["terminate"] = ""

    # let's check whether the trajectory is worth optimising by checking confidence, clashes, and contacts
    # check clashes
    #clash_interface = calculate_clash_score(model_pdb_path, 2.4)
    ca_clashes = calculate_clash_score(model_pdb_path, 2.5, only_ca=True)

    #if clash_interface > 25 or ca_clashes > 0:
    if ca_clashes > 0:
        af_model.aux["log"]["terminate"] = "Clashing"
        update_failures(failure_csv, 'Trajectory_Clashes')
        print("Severe clashes detected, skipping analysis and MPNN optimisation")
        print("")
    else:
        # check if low quality prediction
        if final_plddt < 0.7:
            af_model.aux["log"]["terminate"] = "LowConfidence"
            update_failures(failure_csv, 'Trajectory_final_pLDDT')
            print("Trajectory starting confidence low, skipping analysis and MPNN optimisation")
            print("")
        else:
            # does it have enough contacts to consider?
            binder_contacts = hotspot_residues(model_pdb_path)
            binder_contacts_n = len(binder_contacts.items())

            # if less than 3 contacts then protein is floating above and is not binder
            if binder_contacts_n < 3:
                af_model.aux["log"]["terminate"] = "LowConfidence"
                update_failures(failure_csv, 'Trajectory_Contacts')
                print("Too few contacts at the interface, skipping analysis and MPNN optimisation")
                print("")
            else:
                # phew, trajectory is okay! We can continue
                af_model.aux["log"]["terminate"] = ""
                print("Trajectory successful, final pLDDT: "+str(final_plddt))

    # move low quality prediction:
    if af_model.aux["log"]["terminate"] != "":
        shutil.move(model_pdb_path, design_paths[f"Trajectory/{af_model.aux['log']['terminate']}"])

    ### get the sampled sequence for plotting
    af_model.get_seqs()
    if advanced_settings.get("save_design_trajectory_plots", True):
        plot_trajectory(af_model, design_name, design_paths)

    ### save the hallucination trajectory animation
    if advanced_settings.get("save_design_animations", True):
        save_design_animation(af_model, os.path.join(
            design_paths["Trajectory/Animation"], design_name+".html"))

    if advanced_settings["save_trajectory_pickle"]:
        with open(os.path.join(design_paths["Trajectory/Pickle"], design_name+".pickle"), 'wb') as handle:
            pickle.dump(af_model.aux['all'], handle, protocol=pickle.HIGHEST_PROTOCOL)

    # Calculate and display ipSAE for completed trajectory
    if af_model.aux["log"]["terminate"] == "":
        try:
            traj_pae_matrix = np.array(af_model._tmp["best"]["aux"]["pae"])
            traj_ipsae = calculate_ipsae(traj_pae_matrix, af_model._target_len, length)
            af_model.aux["log"]["ipSAE"] = traj_ipsae['ipSAE']
            print(f"Trajectory ipSAE: {traj_ipsae['ipSAE']:.4f}")
        except Exception as e:
            af_model.aux["log"]["ipSAE"] = None
            vprint(f"[AF2] Could not calculate trajectory ipSAE: {e}")
    else:
        af_model.aux["log"]["ipSAE"] = None

    return af_model

# run prediction for binder with masked template target
def predict_binder_complex(prediction_model, binder_sequence, mpnn_design_name, target_pdb, chain, length, trajectory_pdb, prediction_models, advanced_settings, filters, design_paths, failure_csv, seed=None, use_pyrosetta=True):
    prediction_stats = {}

    # clean sequence
    binder_sequence = re.sub("[^A-Z]", "", binder_sequence.upper())

    # reset filtering conditionals
    pass_af2_filters = True
    filter_failures = {}

    # Removed cyclize_peptide offset in v1.5.2

    # start prediction per AF2 model, 2 are used by default due to masked templates
    for model_num in prediction_models:
        # check to make sure prediction does not exist already
        complex_pdb = os.path.join(design_paths["MPNN"], f"{mpnn_design_name}_model{model_num+1}.pdb")
        if not os.path.exists(complex_pdb):
            # predict model
            prediction_model.predict(seq=binder_sequence, models=[model_num], num_recycles=advanced_settings["num_recycles_validation"], verbose=False)
            prediction_model.save_pdb(complex_pdb)
            prediction_metrics = copy_dict(prediction_model.aux["log"]) # contains plddt, ptm, i_ptm, pae, i_pae

            # Calculate ipSAE from full PAE matrix
            try:
                pae_matrix = np.array(prediction_model.aux["pae"])
                target_len = prediction_model._target_len
                binder_len = prediction_model._binder_len
                ipsae_result = calculate_ipsae(pae_matrix, target_len, binder_len)
                ipsae_score = ipsae_result['ipSAE']
            except Exception:
                ipsae_score = None

            # extract the statistics for the model
            stats = {
                'pLDDT': round(prediction_metrics['plddt'], 2), 
                'pTM': round(prediction_metrics['ptm'], 2), 
                'i_pTM': round(prediction_metrics['i_ptm'], 2), 
                'pAE': round(prediction_metrics['pae'], 2), 
                'i_pAE': round(prediction_metrics['i_pae'], 2),
                'ipSAE': ipsae_score
            }
            prediction_stats[model_num+1] = stats

            # List of filter conditions and corresponding keys
            filter_conditions = [
                (f"{model_num+1}_pLDDT", 'plddt', '>='),
                (f"{model_num+1}_pTM", 'ptm', '>='),
                (f"{model_num+1}_i_pTM", 'i_ptm', '>='),
                (f"{model_num+1}_pAE", 'pae', '<='),
                (f"{model_num+1}_i_pAE", 'i_pae', '<='),
            ]

            # perform initial AF2 values filtering to determine whether to skip relaxation and interface scoring
            for filter_name, metric_key, comparison in filter_conditions:
                threshold = filters.get(filter_name, {}).get("threshold")
                if threshold is not None:
                    if comparison == '>=' and prediction_metrics[metric_key] < threshold:
                        pass_af2_filters = False
                        filter_failures[filter_name] = filter_failures.get(filter_name, 0) + 1
                    elif comparison == '<=' and prediction_metrics[metric_key] > threshold:
                        pass_af2_filters = False
                        filter_failures[filter_name] = filter_failures.get(filter_name, 0) + 1

            if not pass_af2_filters:
                break

    # Update the CSV file with the failure counts
    if filter_failures:
        update_failures(failure_csv, filter_failures)

    # AF2 filters passed, contuing with relaxation
    for model_num in prediction_models:
        complex_pdb = os.path.join(design_paths["MPNN"], f"{mpnn_design_name}_model{model_num+1}.pdb")
        if pass_af2_filters:
            mpnn_relaxed = os.path.join(design_paths["MPNN/Relaxed"], f"{mpnn_design_name}_model{model_num+1}.pdb")
            pr_relax(complex_pdb, mpnn_relaxed, use_pyrosetta=use_pyrosetta)
        else:
            if os.path.exists(complex_pdb):
                os.remove(complex_pdb)

    return prediction_stats, pass_af2_filters, filter_failures

# run prediction for binder alone
def predict_binder_alone(prediction_model, binder_sequence, mpnn_design_name, length, trajectory_pdb, binder_chain, prediction_models, advanced_settings, design_paths, seed=None, use_pyrosetta=True):
    binder_stats = {}

    # prepare sequence for prediction
    binder_sequence = re.sub("[^A-Z]", "", binder_sequence.upper())
    prediction_model.set_seq(binder_sequence)

    # Removed cyclize_peptide offset in v1.5.2

    # predict each model separately
    for model_num in prediction_models:
        # check to make sure prediction does not exist already
        binder_alone_pdb = os.path.join(design_paths["MPNN/Binder"], f"{mpnn_design_name}_model{model_num+1}.pdb")
        if not os.path.exists(binder_alone_pdb):
            # predict model
            prediction_model.predict(models=[model_num], num_recycles=advanced_settings["num_recycles_validation"], verbose=False)
            prediction_model.save_pdb(binder_alone_pdb)
            prediction_metrics = copy_dict(prediction_model.aux["log"]) # contains plddt, ptm, pae

            # align binder model to trajectory binder
            align_pdbs(trajectory_pdb, binder_alone_pdb, binder_chain, "A", use_pyrosetta=use_pyrosetta)

            # extract the statistics for the model
            stats = {
                'pLDDT': round(prediction_metrics['plddt'], 2), 
                'pTM': round(prediction_metrics['ptm'], 2), 
                'pAE': round(prediction_metrics['pae'], 2)
            }
            binder_stats[model_num+1] = stats

    return binder_stats

# run MPNN to generate sequences for binders
def mpnn_gen_sequence(trajectory_pdb, binder_chain, trajectory_interface_residues, advanced_settings, *, clear_gpu_memory=True):
    if clear_gpu_memory:
        clear_mem()

    # initialise MPNN model
    mpnn_model = mk_mpnn_model(backbone_noise=advanced_settings["backbone_noise"], model_name=advanced_settings["model_path"], weights=advanced_settings["mpnn_weights"])

    # check whether keep the interface generated by the trajectory or whether to redesign with MPNN
    design_chains = 'A,' + binder_chain

    if advanced_settings["mpnn_fix_interface"]:
        fixed_positions = 'A,' + trajectory_interface_residues
        fixed_positions = fixed_positions.rstrip(",")
        print("Fixing interface residues: "+trajectory_interface_residues)
    else:
        fixed_positions = 'A'

    # prepare inputs for MPNN
    mpnn_model.prep_inputs(pdb_filename=trajectory_pdb, chain=design_chains, fix_pos=fixed_positions, rm_aa=advanced_settings["omit_AAs"])

    # sample MPNN sequences in parallel
    mpnn_sequences = mpnn_model.sample(temperature=advanced_settings["sampling_temp"], num=1, batch=advanced_settings["num_seqs"]) 

    return mpnn_sequences

# Get pLDDT of best model
def get_best_plddt(af_model, length):
    return round(np.mean(af_model._tmp["best"]["aux"]["plddt"][-length:]),2)

# Define radius of gyration loss for colabdesign
def add_rg_loss(self, weight=0.1):
    '''add radius of gyration loss'''
    def loss_fn(inputs, outputs):
        xyz = outputs["structure_module"]
        ca = xyz["final_atom_positions"][:,residue_constants.atom_order["CA"]]
        ca = ca[-self._binder_len:]
        rg = jnp.sqrt(jnp.square(ca - ca.mean(0)).sum(-1).mean() + 1e-8)
        rg_th = 2.38 * ca.shape[0] ** 0.365

        rg = jax.nn.elu(rg - rg_th)
        return {"rg":rg}

    self._callbacks["model"]["loss"].append(loss_fn)
    self.opt["weights"]["rg"] = weight

# Define interface pTM loss for colabdesign
def add_i_ptm_loss(self, weight=0.1):
    def loss_iptm(inputs, outputs):
        p = 1 - get_ptm(inputs, outputs, interface=True)
        i_ptm = mask_loss(p)
        return {"i_ptm": i_ptm}
    
    self._callbacks["model"]["loss"].append(loss_iptm)
    self.opt["weights"]["i_ptm"] = weight

# add helicity loss
def add_helix_loss(self, weight=0):
    def binder_helicity(inputs, outputs):  
      if "offset" in inputs:
        offset = inputs["offset"]
      else:
        idx = inputs["residue_index"].flatten()
        offset = idx[:,None] - idx[None,:]

      # define distogram
      dgram = outputs["distogram"]["logits"]
      dgram_bins = get_dgram_bins(outputs)
      mask_2d = np.outer(np.append(np.zeros(self._target_len), np.ones(self._binder_len)), np.append(np.zeros(self._target_len), np.ones(self._binder_len)))

      x = _get_con_loss(dgram, dgram_bins, cutoff=6.0, binary=True)
      if offset is None:
        if mask_2d is None:
          helix_loss = jnp.diagonal(x,3).mean()
        else:
          helix_loss = jnp.diagonal(x * mask_2d,3).sum() + (jnp.diagonal(mask_2d,3).sum() + 1e-8)
      else:
        mask = offset == 3
        if mask_2d is not None:
          mask = jnp.where(mask_2d,mask,0)
        helix_loss = jnp.where(mask,x,0.0).sum() / (mask.sum() + 1e-8)

      return {"helix":helix_loss}
    self._callbacks["model"]["loss"].append(binder_helicity)
    self.opt["weights"]["helix"] = weight

# add N- and C-terminus distance loss
def add_termini_distance_loss(self, weight=0.1, threshold_distance=7.0):
    '''Add loss penalizing the distance between N and C termini'''
    def loss_fn(inputs, outputs):
        xyz = outputs["structure_module"]
        ca = xyz["final_atom_positions"][:, residue_constants.atom_order["CA"]]
        ca = ca[-self._binder_len:]  # Considering only the last _binder_len residues

        # Extract N-terminus (first CA atom) and C-terminus (last CA atom)
        n_terminus = ca[0]
        c_terminus = ca[-1]

        # Compute the distance between N and C termini
        termini_distance = jnp.linalg.norm(n_terminus - c_terminus)

        # Compute the deviation from the threshold distance using ELU activation
        deviation = jax.nn.elu(termini_distance - threshold_distance)

        # Ensure the loss is never lower than 0
        termini_distance_loss = jax.nn.relu(deviation)
        return {"NC": termini_distance_loss}

    # Append the loss function to the model callbacks
    self._callbacks["model"]["loss"].append(loss_fn)
    self.opt["weights"]["NC"] = weight

# plot design trajectory losses
def plot_trajectory(af_model, design_name, design_paths):
    metrics_to_plot = ['loss', 'plddt', 'ptm', 'i_ptm', 'con', 'i_con', 'pae', 'i_pae', 'rg', 'mpnn']
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

    for index, metric in enumerate(metrics_to_plot):
        if metric in af_model.aux["log"]:
            # Create a new figure for each metric
            plt.figure()

            loss = af_model.get_loss(metric)
            # Create an x axis for iterations
            iterations = range(1, len(loss) + 1)

            plt.plot(iterations, loss, label=f'{metric}', color=colors[index % len(colors)])

            # Add labels and a legend
            plt.xlabel('Iterations')
            plt.ylabel(metric)
            plt.title(design_name)
            plt.legend()
            plt.grid(True)

            # Save the plot
            plt.savefig(os.path.join(design_paths["Trajectory/Plots"], design_name+"_"+metric+".png"), dpi=150)
            
            # Close the figure
            plt.close()
