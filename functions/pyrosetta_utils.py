####################################
################ PyRosetta functions
####################################
### Import dependencies
import os
import time
from .generic_utils import clean_pdb
from .logging_utils import vprint
from .biopython_utils import hotspot_residues, biopython_unaligned_rmsd, biopython_align_pdbs
from . import pr_alternative_utils as alt

# Conditionally import PyRosetta - will be available if initialized successfully
pr = None
try:
    import pyrosetta as pr
    from pyrosetta.rosetta.core.kinematics import MoveMap
    from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
    from pyrosetta.rosetta.protocols.simple_moves import AlignChainMover
    from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.core.simple_metrics.metrics import RMSDMetric
    from pyrosetta.rosetta.core.select import get_residues_from_subset
    from pyrosetta.rosetta.core.io import pose_from_pose
    from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
    PYROSETTA_AVAILABLE = True
except ImportError:
    PYROSETTA_AVAILABLE = False
    # Suppress import-time warnings; runtime messaging is handled in bindcraft.py


def score_interface(pdb_file, binder_chain="B", use_pyrosetta=True):
    """
    Calculate interface scores using PyRosetta or alternative methods.
    Dispatches to appropriate implementation based on availability.
    """
    basename = os.path.basename(pdb_file)
    t0_all = time.time()
    try:
        # Handle PyRosetta-free mode
        if not use_pyrosetta or not PYROSETTA_AVAILABLE:
            return alt.pr_alternative_score_interface(pdb_file, binder_chain)
    except ImportError as e:
        if "pr_alternative_utils" in str(e):
            raise ImportError("Failed to import alternative implementations module") from e
        raise
        
    # Regular PyRosetta mode
    vprint(f"[Rosetta-Score] Initiating PyRosetta scoring for {basename} (binder={binder_chain})")
    # load pose
    t0_pose = time.time()
    pose = pr.pose_from_pdb(pdb_file)
    vprint(f"[Rosetta-Score] Loaded pose in {time.time()-t0_pose:.2f}s")

    # analyze interface statistics
    iam = InterfaceAnalyzerMover()
    iam.set_interface("A_B")
    scorefxn = pr.get_fa_scorefxn()
    iam.set_scorefunction(scorefxn)
    iam.set_compute_packstat(True)
    iam.set_compute_interface_energy(True)
    iam.set_calc_dSASA(True)
    iam.set_calc_hbond_sasaE(True)
    iam.set_compute_interface_sc(True)
    iam.set_pack_separated(True)
    t0_iam = time.time()
    iam.apply(pose)
    vprint(f"[Rosetta-Score] InterfaceAnalyzerMover applied in {time.time()-t0_iam:.2f}s")

    # Initialize dictionary with all amino acids
    interface_AA = {aa: 0 for aa in 'ACDEFGHIKLMNPQRSTVWY'}

    # Initialize list to store PDB residue IDs at the interface
    t0_if = time.time()
    interface_residues_set = hotspot_residues(pdb_file, binder_chain)
    interface_residues_pdb_ids = []
    
    # Iterate over the interface residues
    for pdb_res_num, aa_type in interface_residues_set.items():
        # Increase the count for this amino acid type
        interface_AA[aa_type] += 1

        # Append the binder_chain and the PDB residue number to the list
        interface_residues_pdb_ids.append(f"{binder_chain}{pdb_res_num}")
    vprint(f"[Rosetta-Score] Found {len(interface_residues_pdb_ids)} interface residues in {time.time()-t0_if:.2f}s")

    # count interface residues
    interface_nres = len(interface_residues_pdb_ids)

    # Convert the list into a comma-separated string
    interface_residues_pdb_ids_str = ','.join(interface_residues_pdb_ids)

    # Calculate the percentage of hydrophobic residues at the interface of the binder
    hydrophobic_aa = set('ACFGILMPVWY')
    hydrophobic_count = sum(interface_AA[aa] for aa in hydrophobic_aa)
    if interface_nres != 0:
        interface_hydrophobicity = (hydrophobic_count / interface_nres) * 100
    else:
        interface_hydrophobicity = 0

    # retrieve statistics
    interfacescore = iam.get_all_data()
    interface_sc = interfacescore.sc_value # shape complementarity
    interface_interface_hbonds = interfacescore.interface_hbonds # number of interface H-bonds
    interface_dG = iam.get_interface_dG() # interface dG
    interface_dSASA = iam.get_interface_delta_sasa() # interface dSASA (interface surface area)
    interface_packstat = iam.get_interface_packstat() # interface pack stat score
    interface_dG_SASA_ratio = interfacescore.dG_dSASA_ratio * 100 # ratio of dG/dSASA (normalised energy for interface area size)
    buns_filter = XmlObjects.static_get_filter('<BuriedUnsatHbonds report_all_heavy_atom_unsats="true" scorefxn="scorefxn" ignore_surface_res="false" use_ddG_style="true" dalphaball_sasa="1" probe_radius="1.1" burial_cutoff_apo="0.2" confidence="0" />')
    interface_delta_unsat_hbonds = buns_filter.report_sm(pose)

    if interface_nres != 0:
        interface_hbond_percentage = (interface_interface_hbonds / interface_nres) * 100 # Hbonds per interface size percentage
        interface_bunsch_percentage = (interface_delta_unsat_hbonds / interface_nres) * 100 # Unsaturated H-bonds per percentage
    else:
        interface_hbond_percentage = None
        interface_bunsch_percentage = None

    # calculate binder energy score
    chain_design = ChainSelector(binder_chain)
    t0_energy = time.time()
    tem = pr.rosetta.core.simple_metrics.metrics.TotalEnergyMetric()
    tem.set_scorefunction(scorefxn)
    tem.set_residue_selector(chain_design)
    binder_score = tem.calculate(pose)
    vprint(f"[Rosetta-Score] TotalEnergyMetric calculated in {time.time()-t0_energy:.2f}s")

    # calculate binder SASA fraction
    t0_sasa = time.time()
    bsasa = pr.rosetta.core.simple_metrics.metrics.SasaMetric()
    bsasa.set_residue_selector(chain_design)
    binder_sasa = bsasa.calculate(pose)
    vprint(f"[Rosetta-Score] SasaMetric calculated in {time.time()-t0_sasa:.2f}s")

    if binder_sasa > 0:
        interface_binder_fraction = (interface_dSASA / binder_sasa) * 100
    else:
        interface_binder_fraction = 0

    # calculate surface hydrophobicity
    binder_pose = {pose.pdb_info().chain(pose.conformation().chain_begin(i)): p for i, p in zip(range(1, pose.num_chains()+1), pose.split_by_chain())}[binder_chain]

    t0_layer = time.time()
    layer_sel = pr.rosetta.core.select.residue_selector.LayerSelector()
    layer_sel.set_layers(pick_core = False, pick_boundary = False, pick_surface = True)
    surface_res = layer_sel.apply(binder_pose)
    vprint(f"[Rosetta-Score] Surface layer selection in {time.time()-t0_layer:.2f}s")

    exp_apol_count = 0
    total_count = 0 
    
    # count apolar and aromatic residues at the surface
    for i in range(1, len(surface_res) + 1):
        if surface_res[i]:
            res = binder_pose.residue(i)

            # count apolar and aromatic residues as hydrophobic
            if res.is_apolar() or res.name() == 'PHE' or res.name() == 'TRP' or res.name() == 'TYR':
                exp_apol_count += 1
            total_count += 1

    surface_hydrophobicity = exp_apol_count/total_count if total_count > 0 else 0.0 # Added safety for division by zero

    # output interface score array and amino acid counts at the interface
    interface_scores = {
    'binder_score': binder_score,
    'surface_hydrophobicity': surface_hydrophobicity,
    'interface_sc': interface_sc,
    'interface_packstat': interface_packstat,
    'interface_dG': interface_dG,
    'interface_dSASA': interface_dSASA,
    'interface_dG_SASA_ratio': interface_dG_SASA_ratio,
    'interface_fraction': interface_binder_fraction,
    'interface_hydrophobicity': interface_hydrophobicity,
    'interface_nres': interface_nres,
    'interface_interface_hbonds': interface_interface_hbonds,
    'interface_hbond_percentage': interface_hbond_percentage,
    'interface_delta_unsat_hbonds': interface_delta_unsat_hbonds,
    'interface_delta_unsat_hbonds_percentage': interface_bunsch_percentage
    }

    # round to two decimal places
    interface_scores = {k: round(v, 2) if isinstance(v, float) else v for k, v in interface_scores.items()}

    vprint(f"[Rosetta-Score] Completed PyRosetta scoring for {basename} in {time.time()-t0_all:.2f}s")
    return interface_scores, interface_AA, interface_residues_pdb_ids_str

def align_pdbs(reference_pdb, align_pdb, reference_chain_id, align_chain_id, use_pyrosetta=True):
    """
    Align PDB structures using PyRosetta or Biopython.
    """
    if not use_pyrosetta or not PYROSETTA_AVAILABLE:
        biopython_align_pdbs(reference_pdb, align_pdb, reference_chain_id, align_chain_id)
        return
        
    # initiate poses
    t0 = time.time()
    basename = os.path.basename(align_pdb)
    vprint(f"[Rosetta-Align] Initiating alignment for {basename} (ref_chain={reference_chain_id}, align_chain={align_chain_id})")
    reference_pose = pr.pose_from_pdb(reference_pdb)
    align_pose = pr.pose_from_pdb(align_pdb)

    align = AlignChainMover()
    align.pose(reference_pose)

    # If the chain IDs contain commas, split them and only take the first value
    reference_chain_id = reference_chain_id.split(',')[0]
    align_chain_id = align_chain_id.split(',')[0]

    # Get the chain number corresponding to the chain ID in the poses
    reference_chain = pr.rosetta.core.pose.get_chain_id_from_chain(reference_chain_id, reference_pose)
    align_chain = pr.rosetta.core.pose.get_chain_id_from_chain(align_chain_id, align_pose)

    align.source_chain(align_chain)
    align.target_chain(reference_chain)
    align.apply(align_pose)
    vprint(f"[Rosetta-Align] Applied alignment in {time.time()-t0:.2f}s")

    # Overwrite aligned pdb
    t1 = time.time()
    align_pose.dump_pdb(align_pdb)
    clean_pdb(align_pdb)
    vprint(f"[Rosetta-Align] Saved and cleaned aligned PDB in {time.time()-t1:.2f}s")

def unaligned_rmsd(reference_pdb, align_pdb, reference_chain_id, align_chain_id, use_pyrosetta=True):
    """
    Calculate RMSD without alignment using PyRosetta or Biopython.
    """
    if not use_pyrosetta or not PYROSETTA_AVAILABLE:
        return biopython_unaligned_rmsd(reference_pdb, align_pdb, reference_chain_id, align_chain_id)
        
    t0 = time.time()
    basename = os.path.basename(align_pdb)
    vprint(f"[Rosetta-RMSD] Initiating unaligned RMSD for {basename} (ref_chain={reference_chain_id}, align_chain={align_chain_id})")
    reference_pose = pr.pose_from_pdb(reference_pdb)
    align_pose = pr.pose_from_pdb(align_pdb)

    # Define chain selectors for the reference and align chains
    reference_chain_selector = ChainSelector(reference_chain_id)
    align_chain_selector = ChainSelector(align_chain_id)

    # Apply selectors to get residue subsets
    reference_chain_subset = reference_chain_selector.apply(reference_pose)
    align_chain_subset = align_chain_selector.apply(align_pose)

    # Convert subsets to residue index vectors
    reference_residue_indices = get_residues_from_subset(reference_chain_subset)
    align_residue_indices = get_residues_from_subset(align_chain_subset)

    # Create empty subposes
    reference_chain_pose = pr.Pose()
    align_chain_pose = pr.Pose()

    # Fill subposes
    pose_from_pose(reference_chain_pose, reference_pose, reference_residue_indices)
    pose_from_pose(align_chain_pose, align_pose, align_residue_indices)

    # Calculate RMSD using the RMSDMetric
    rmsd_metric = RMSDMetric()
    rmsd_metric.set_comparison_pose(reference_chain_pose)
    rmsd = rmsd_metric.calculate(align_chain_pose)
    vprint(f"[Rosetta-RMSD] Computed RMSD={rmsd:.2f} in {time.time()-t0:.2f}s")

    return round(rmsd, 2)

def pr_relax(pdb_file, relaxed_pdb_path, use_pyrosetta=True):
    """
    Relax PDB structure using PyRosetta FastRelax or OpenMM.
    """
    if use_pyrosetta and PYROSETTA_AVAILABLE:
        basename = os.path.basename(pdb_file)
        vprint(f"[Rosetta-Relax] Initiating FastRelax for {basename}")
        t0_all = time.time()
        pose = pr.pose_from_pdb(pdb_file)
        start_pose = pose.clone()

        ### Generate movemaps
        mmf = MoveMap()
        mmf.set_chi(True) # enable sidechain movement
        mmf.set_bb(True) # enable backbone movement, can be disabled to increase speed by 30% but makes metrics look worse on average
        mmf.set_jump(False) # disable whole chain movement

        # Run FastRelax
        fastrelax = FastRelax()
        scorefxn = pr.get_fa_scorefxn()
        fastrelax.set_scorefxn(scorefxn)
        fastrelax.set_movemap(mmf) # set MoveMap
        fastrelax.max_iter(200) # default iterations is 2500
        fastrelax.min_type("lbfgs_armijo_nonmonotone")
        fastrelax.constrain_relax_to_start_coords(True)
        t0_relax = time.time()
        fastrelax.apply(pose)
        vprint(f"[Rosetta-Relax] FastRelax applied in {time.time()-t0_relax:.2f}s")

        # Align relaxed structure to original trajectory
        align = AlignChainMover()
        align.source_chain(0)
        align.target_chain(0)
        align.pose(start_pose)
        t0_align = time.time()
        align.apply(pose)
        vprint(f"[Rosetta-Relax] Post-relax alignment applied in {time.time()-t0_align:.2f}s")

        # Copy B factors from start_pose to pose
        t0_bfac = time.time()
        for resid in range(1, pose.total_residue() + 1):
            if pose.residue(resid).is_protein():
                # Get the B factor of the first heavy atom in the residue
                bfactor = start_pose.pdb_info().bfactor(resid, 1)
                for atom_id in range(1, pose.residue(resid).natoms() + 1):
                    pose.pdb_info().bfactor(resid, atom_id, bfactor)
        vprint(f"[Rosetta-Relax] B-factor transfer in {time.time()-t0_bfac:.2f}s")

        # output relaxed and aligned PDB
        t0_save = time.time()
        pose.dump_pdb(relaxed_pdb_path)
        clean_pdb(relaxed_pdb_path)
        vprint(f"[Rosetta-Relax] Saved and cleaned relaxed PDB in {time.time()-t0_save:.2f}s")
        vprint(f"[Rosetta-Relax] Completed FastRelax for {basename} in {time.time()-t0_all:.2f}s")
    else:
        # Quiet: direct OpenMM fallback without extra banner
        openmm_gpu = True # Default to True for GPU usage in OpenMM fallback
        # Run OpenMM relax in a fresh subprocess to fully reset OpenCL context per run
        alt.openmm_relax_subprocess(pdb_file, relaxed_pdb_path, use_gpu_relax=openmm_gpu)
        