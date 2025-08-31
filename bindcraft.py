####################################
###################### BindCraft Run
####################################
### Import dependencies
import gc
from functions import *
from functions.generic_utils import insert_data # Explicit import for insert_data
from functions.biopython_utils import clear_dssp_cache # Explicit import for DSSP cache management
import logging
import os
import sys
import subprocess
try:
    import resource  # POSIX-only; used to raise RLIMIT_NOFILE (ulimit -n)
except Exception:
    resource = None

def _bump_open_files_limit(min_soft=65536):
    """Attempt to raise the soft RLIMIT_NOFILE up to min_soft (not above hard)."""
    if resource is None:
        print("Warning: 'resource' module not available; cannot adjust open files limit.")
        return
    try:
        soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        desired_soft = min(max(soft, int(min_soft)), hard if hard != resource.RLIM_INFINITY else max(soft, int(min_soft)))
        if desired_soft > soft:
            resource.setrlimit(resource.RLIMIT_NOFILE, (desired_soft, hard))
            print(f"Adjusted open files soft limit: {soft} -> {desired_soft} (hard={hard})")
        else:
            print(f"Open files limits OK (soft={soft}, hard={hard})")
    except Exception as e:
        print(f"Warning: Unable to adjust open files limit: {e}")

# Raise file descriptor soft limit early to avoid 'Too many open files'
_bump_open_files_limit(min_soft=65536)

# Defer GPU availability check until after CLI/interactive handling

def ensure_binaries_executable(use_pyrosetta=True):
    """Ensure all required binaries in functions/ are executable."""
    bindcraft_folder = os.path.dirname(os.path.realpath(__file__))
    functions_dir = os.path.join(bindcraft_folder, "functions")
    
    # Always needed binaries
    binaries = ["dssp", "sc"]
    
    # Only add DAlphaBall.gcc if using PyRosetta
    if use_pyrosetta:
        binaries.append("DAlphaBall.gcc")
    
    for binary in binaries:
        binary_path = os.path.join(functions_dir, binary)
        if os.path.isfile(binary_path):
            try:
                # Check if already executable
                if not os.access(binary_path, os.X_OK):
                    # Make executable
                    current_mode = os.stat(binary_path).st_mode
                    os.chmod(binary_path, current_mode | 0o755)
                    print(f"Made {binary} executable")
            except Exception as e:
                print(f"Warning: Failed to make {binary} executable: {e}")

# Ensure binaries are executable at startup (will be called again with proper use_pyrosetta flag later)
ensure_binaries_executable()

######################################
### parse input paths and interactive mode
parser = argparse.ArgumentParser(description='Script to run BindCraft binder design.')

parser.add_argument('--settings', '-s', type=str, required=False,
                    help='Path to the basic settings.json file. If omitted in a TTY, interactive mode is used.')
parser.add_argument('--filters', '-f', type=str, default='./settings_filters/default_filters.json',
                    help='Path to the filters.json file used to filter design. If not provided, default will be used.')
parser.add_argument('--advanced', '-a', type=str, default='./settings_advanced/default_4stage_multimer.json',
                    help='Path to the advanced.json file with additional design settings. If not provided, default will be used.')
parser.add_argument('--no-pyrosetta', action='store_true',
                    help='Run without PyRosetta (skips relaxation and PyRosetta-based scoring)')
parser.add_argument('--verbose', action='store_true',
                    help='Enable detailed timing/progress logs')
parser.add_argument('--no-plots', action='store_true',
                    help='Disable saving design trajectory plots (overrides advanced settings)')
parser.add_argument('--no-animations', action='store_true',
                    help='Disable saving design animations (overrides advanced settings)')
parser.add_argument('--interactive', action='store_true',
                    help='Force interactive mode to collect target settings and options')

args = parser.parse_args()

def _isatty_stdin():
    try:
        return sys.stdin.isatty()
    except Exception:
        return False

def _input_with_default(prompt_text, default_value=None):
    if default_value is None:
        return input(prompt_text).strip()
    resp = input(f"{prompt_text} ").strip()
    return resp if resp else default_value

def _yes_no(prompt_text, default_yes=False):
    default_hint = 'Y/n' if default_yes else 'y/N'
    resp = input(f"{prompt_text} ({default_hint}): ").strip().lower()
    if resp == '':
        return default_yes
    return resp in ('y', 'yes')

def _list_json_choices(folder_path):
    try:
        entries = [f for f in os.listdir(folder_path) if f.endswith('.json') and not f.startswith('.')]
    except Exception:
        entries = []
    entries.sort()
    # return list of (display_name_without_ext, abspath)
    return [(os.path.splitext(f)[0], os.path.join(folder_path, f)) for f in entries]

def _prompt_interactive_and_prepare_args(args):
    base_dir = os.path.dirname(os.path.realpath(__file__))
    filters_dir = os.path.join(base_dir, 'settings_filters')
    advanced_dir = os.path.join(base_dir, 'settings_advanced')

    print("\nBindCraft Interactive Setup\n")
    # Container plan state shared across confirmation
    container_plan = {
        "use": False,
        "image": None,
        "gpu_idx": None
    }

    while True:
        # Design type selection
        print("Design type:")
        print("1. Miniprotein (31+ aa)")
        print("2. Peptide (8-30 aa)")
        dtype_choice = _input_with_default("Choose design type (press Enter for Miniprotein):", "")
        is_peptide = (dtype_choice.strip() == '2')

        # Required inputs
        project_name = _input_with_default("Enter project/binder name:", None)
        while not project_name:
            print("Project name is required.")
            project_name = _input_with_default("Enter project/binder name:", None)

        # Require a valid existing PDB file path; re-prompt until valid
        while True:
            pdb_raw = _input_with_default("Enter path to PDB file:", None)
            if not pdb_raw:
                print("PDB path is required.")
                continue
            candidate = os.path.abspath(os.path.expanduser(pdb_raw))
            if os.path.isfile(candidate):
                pdb_path = candidate
                break
            print(f"Error: No PDB file found at '{candidate}'. Please re-enter.")

        output_dir = _input_with_default("Enter output directory:", os.path.join(os.getcwd(), f"{project_name}_bindcraft_out"))
        output_dir = os.path.abspath(os.path.expanduser(output_dir))
        os.makedirs(output_dir, exist_ok=True)

        chains = _input_with_default("Enter target chains (e.g., A or A,B):", "A")
        hotspots = _input_with_default("Enter hotspot residue(s) for BindCraft to target. Use format: chain letter + residue numbers (e.g., 'A1,B20-25'). Leave empty for no preference:", "")

        if is_peptide:
            lengths_prompt_default = "8 25"
            lengths_input = _input_with_default("Enter peptide min and max lengths (8-30) separated by space or comma (default 8 25):", lengths_prompt_default)
        else:
            lengths_prompt_default = "65 150"
            lengths_input = _input_with_default("Enter miniprotein min and max lengths separated by space or comma (min>=31, default 65 150):", lengths_prompt_default)
        try:
            normalized = lengths_input.replace(',', ' ').split()
            min_len_val, max_len_val = int(normalized[0]), int(normalized[1])
        except Exception:
            if is_peptide:
                min_len_val, max_len_val = 8, 25
            else:
                min_len_val, max_len_val = 65, 150
        if is_peptide:
            # Clamp within [8,30]
            min_len_val = max(8, min(min_len_val, 30))
            max_len_val = max(8, min(max_len_val, 30))
            if min_len_val > max_len_val:
                min_len_val, max_len_val = max_len_val, min_len_val
        else:
            # Enforce min >= 31; ensure order
            if min_len_val < 31:
                min_len_val = 31
            if max_len_val < min_len_val:
                max_len_val = min_len_val
        lengths = [min_len_val, max_len_val]

        num_designs_str = _input_with_default("Enter number of final designs (default 100):", "100")
        try:
            num_designs = int(num_designs_str)
        except Exception:
            num_designs = 100

        # List choices
        print("\nAvailable filter settings:")
        filter_choices_all = _list_json_choices(filters_dir)
        name_to_filter = {name: path for name, path in filter_choices_all}
        if is_peptide:
            filter_order = ['peptide_filters', 'peptide_relaxed_filters', 'no_filters']
            default_filter_name = 'peptide_filters'
        else:
            filter_order = ['default_filters', 'relaxed_filters', 'no_filters']
            default_filter_name = 'default_filters'
        ordered_filters = [(name, name_to_filter[name]) for name in filter_order if name in name_to_filter]
        for i, (name, _) in enumerate(ordered_filters, 1):
            print(f"{i}. {name}")
        filter_idx = _input_with_default(f"Choose filter (press Enter for {default_filter_name}):", "")
        if filter_idx:
            try:
                filter_idx_int = int(filter_idx)
                selected_filter = ordered_filters[filter_idx_int - 1][1]
            except Exception:
                selected_filter = name_to_filter.get(default_filter_name, os.path.join(filters_dir, f"{default_filter_name}.json"))
        else:
            selected_filter = name_to_filter.get(default_filter_name, os.path.join(filters_dir, f"{default_filter_name}.json"))

        print("\nAvailable advanced settings:")
        advanced_choices_all = _list_json_choices(advanced_dir)
        name_to_adv = {name: path for name, path in advanced_choices_all}
        if is_peptide:
            adv_order = [
                'peptide_3stage_multimer',
                'peptide_3stage_multimer_mpnn',
                'peptide_3stage_multimer_flexible',
                'peptide_3stage_multimer_mpnn_flexible'
            ]
            default_adv_name = 'peptide_3stage_multimer'
        else:
            adv_order = [
                'default_4stage_multimer',
                'default_4stage_multimer_mpnn',
                'default_4stage_multimer_flexible',
                'default_4stage_multimer_hardtarget',
                'default_4stage_multimer_flexible_hardtarget',
                'default_4stage_multimer_mpnn_flexible',
                'default_4stage_multimer_mpnn_hardtarget',
                'default_4stage_multimer_mpnn_flexible_hardtarget',
                'betasheet_4stage_multimer',
                'betasheet_4stage_multimer_mpnn',
                'betasheet_4stage_multimer_flexible',
                'betasheet_4stage_multimer_hardtarget',
                'betasheet_4stage_multimer_flexible_hardtarget',
                'betasheet_4stage_multimer_mpnn_flexible',
                'betasheet_4stage_multimer_mpnn_hardtarget',
                'betasheet_4stage_multimer_mpnn_flexible_hardtarget'
            ]
            default_adv_name = 'default_4stage_multimer'
        ordered_adv = [(name, name_to_adv[name]) for name in adv_order if name in name_to_adv]
        for i, (name, _) in enumerate(ordered_adv, 1):
            print(f"{i}. {name}")
        advanced_idx = _input_with_default(f"Choose advanced (press Enter for {default_adv_name}):", "")
        if advanced_idx:
            try:
                advanced_idx_int = int(advanced_idx)
                selected_advanced = ordered_adv[advanced_idx_int - 1][1]
            except Exception:
                selected_advanced = name_to_adv.get(default_adv_name, os.path.join(advanced_dir, f"{default_adv_name}.json"))
        else:
            selected_advanced = name_to_adv.get(default_adv_name, os.path.join(advanced_dir, f"{default_adv_name}.json"))

        # Toggles
        verbose = _yes_no("Enable verbose output?", default_yes=False)
        plots_on = _yes_no("Enable saving plots?", default_yes=True)
        animations_on = _yes_no("Enable saving animations?", default_yes=True)
        run_with_pyrosetta = _yes_no("Run with PyRosetta?", default_yes=True)

        # Container selection BEFORE final confirmation
        container_plan = {k: None for k in container_plan}  # reset plan
        container_plan["use"] = _yes_no("Run using a Docker container?", default_yes=False)
        if container_plan["use"]:
            # Launch new container from image only
            container_plan["image"] = _input_with_default("Docker image to use (default bindcraft:latest):", "bindcraft:latest")
            container_plan["gpu_idx"] = _input_with_default("GPU index to expose (default 0):", "0")

        # Summary for confirmation
        print("\nConfiguration Summary:")
        print(f"Project Name: {project_name}")
        print(f"PDB File: {pdb_path}")
        print(f"Output Directory: {output_dir}")
        print(f"Chains: {chains}")
        print(f"Hotspots: {hotspots if hotspots else 'None'}")
        print(f"Length Range: {lengths}")
        print(f"Design Type: {'Peptide' if is_peptide else 'Miniprotein'}")
        print(f"Number of Final Designs: {num_designs}")
        print(f"Filter Setting: {os.path.splitext(os.path.basename(selected_filter))[0]}")
        print(f"Advanced Setting: {os.path.splitext(os.path.basename(selected_advanced))[0]}")
        print(f"Verbose: {'Yes' if verbose else 'No'}")
        print(f"Plots: {'On' if plots_on else 'Off'}")
        print(f"Animations: {'On' if animations_on else 'Off'}")
        print(f"PyRosetta: {'On' if run_with_pyrosetta else 'Off'}")
        if container_plan["use"]:
            print(f"Docker: Yes (launch new container) -> image={container_plan['image']} gpu={container_plan['gpu_idx']}")
        else:
            print("Docker: No")

        if _yes_no("Proceed with these settings?", default_yes=True):
            break
        else:
            print("Let's re-enter the details.\n")

    # Prepare target settings JSON
    target_settings = {
        "design_path": output_dir,
        "binder_name": project_name,
        "starting_pdb": pdb_path,
        "chains": chains,
        "target_hotspot_residues": hotspots,
        "lengths": lengths,
        "number_of_final_designs": num_designs
    }

    timestamp = time.strftime('%Y%m%d_%H%M%S')
    settings_filename = f"interactive_{project_name}_{timestamp}.json"
    settings_path_out = os.path.join(output_dir, settings_filename)
    try:
        with open(settings_path_out, 'w') as f:
            json.dump(target_settings, f, indent=2)
    except Exception as e:
        print(f"Error writing settings JSON: {e}")
        sys.exit(1)

    # Map to args
    args.settings = settings_path_out
    args.filters = selected_filter
    args.advanced = selected_advanced
    args.verbose = verbose
    args.no_plots = (not plots_on)
    args.no_animations = (not animations_on)
    args.no_pyrosetta = (not run_with_pyrosetta)

    # Execute container plan if selected
    if container_plan["use"] and container_plan["image"]:
        docker_image = container_plan["image"]
        gpu_idx = container_plan["gpu_idx"] or "0"

        mounts = []
        cwd = os.getcwd()
        mounts.append((cwd, cwd))
        pdb_parent = os.path.dirname(pdb_path)
        if pdb_parent and os.path.abspath(pdb_parent) != cwd:
            mounts.append((pdb_parent, pdb_parent))
        mounts.append((output_dir, output_dir))

        docker_cmd = [
            "docker", "run", "--rm", "-it",
            "--gpus", f"device={gpu_idx}",
        ]
        for host, cont in mounts:
            docker_cmd.extend(["-v", f"{host}:{cont}"])
        docker_cmd.extend(["-w", cwd, docker_image, "python", "bindcraft.py",
                           "-s", args.settings, "-f", args.filters, "-a", args.advanced])
        if args.no_pyrosetta:
            docker_cmd.append("--no-pyrosetta")
        if args.verbose:
            docker_cmd.append("--verbose")
        if args.no_plots:
            docker_cmd.append("--no-plots")
        if args.no_animations:
            docker_cmd.append("--no-animations")
        subprocess.run(docker_cmd, check=False)
        sys.exit(0)

    return args

# Enter interactive mode if requested or if no settings were provided in a TTY
if args.interactive or (not args.settings and _isatty_stdin()):
    args = _prompt_interactive_and_prepare_args(args)
elif not args.settings and not _isatty_stdin():
    # No TTY and no settings -> cannot prompt
    print("Error: --settings is required in non-interactive environments.")
    sys.exit(1)

# perform checks of input setting files
settings_path, filters_path, advanced_path = perform_input_check(args)

# Configure standard logging based on --verbose
logging.basicConfig(
    level=logging.WARNING,
    format='%(asctime)s %(levelname)s %(name)s: %(message)s'
)

# Reduce noise from third-party libraries regardless of verbosity
for noisy_logger in (
    "jax", "jaxlib", "jax._src", "absl", "flax", "colabdesign", "tensorflow", "xla"):
    logging.getLogger(noisy_logger).setLevel(logging.WARNING)

# Enable detailed logs only for our modules when --verbose is set
if args.verbose:
    logging.getLogger("functions").setLevel(logging.DEBUG)
    logging.getLogger("bindcraft").setLevel(logging.DEBUG)

# Check if JAX-capable GPU is available, otherwise exit (after interactive/container handling)
check_jax_gpu()

### load settings from JSON
target_settings, advanced_settings, filters = load_json_settings(settings_path, filters_path, advanced_path)

settings_file = os.path.basename(settings_path).split('.')[0]
filters_file = os.path.basename(filters_path).split('.')[0]
advanced_file = os.path.basename(advanced_path).split('.')[0]

### load AF2 model settings
design_models, prediction_models, multimer_validation = load_af2_models(advanced_settings["use_multimer_design"])

### perform checks on advanced_settings
bindcraft_folder = os.path.dirname(os.path.realpath(__file__))
advanced_settings = perform_advanced_settings_check(advanced_settings, bindcraft_folder)

# CLI overrides for plots/animations
if args.no_plots:
    advanced_settings["save_design_trajectory_plots"] = False
if args.no_animations:
    advanced_settings["save_design_animations"] = False

### generate directories, design path names can be found within the function
design_paths = generate_directories(target_settings["design_path"])

### generate dataframes
trajectory_labels, design_labels, final_labels = generate_dataframe_labels()

trajectory_csv = os.path.join(target_settings["design_path"], 'trajectory_stats.csv')
mpnn_csv = os.path.join(target_settings["design_path"], 'mpnn_design_stats.csv')
final_csv = os.path.join(target_settings["design_path"], 'final_design_stats.csv')
failure_csv = os.path.join(target_settings["design_path"], 'failure_csv.csv')

create_dataframe(trajectory_csv, trajectory_labels)
create_dataframe(mpnn_csv, design_labels)
create_dataframe(final_csv, final_labels)
generate_filter_pass_csv(failure_csv, args.filters)

### Define and initialize rejected_mpnn_full_stats.csv
# Ensure failure_csv exists and has headers to read its column structure
if not os.path.exists(failure_csv):
    # This should ideally not happen if generate_filter_pass_csv worked, create an empty one if it's missing.
    temp_failure_df_for_cols = pd.DataFrame()
    print(f"Warning: {failure_csv} was not found after generate_filter_pass_csv. rejected_mpnn_full_stats.csv might have incorrect filter columns.")
else:
    try:
        temp_failure_df_for_cols = pd.read_csv(failure_csv)
    except pd.errors.EmptyDataError:
        # If failure_csv is empty we need to get column names from how generate_filter_pass_csv would create them.
        print(f"Warning: {failure_csv} is empty. rejected_mpnn_full_stats.csv may lack detailed filter columns initially if no filters are defined or an issue occurred.")
        temp_failure_df_for_cols = pd.DataFrame() # Fallback- we don't want BindCraft to crash over this

filter_column_names_for_rejected_log = temp_failure_df_for_cols.columns.tolist()
del temp_failure_df_for_cols # Free memory

rejected_stats_columns = ['Design', 'Sequence'] + filter_column_names_for_rejected_log
rejected_mpnn_full_stats_csv = os.path.join(target_settings["design_path"], 'rejected_mpnn_full_stats.csv')
create_dataframe(rejected_mpnn_full_stats_csv, rejected_stats_columns)
####################################
####################################
####################################
### initialise PyRosetta if not disabled
use_pyrosetta = False

if args.no_pyrosetta:
    # Quiet when user explicitly disables PyRosetta
    print("Running in PyRosetta-free mode as requested by --no-pyrosetta flag.")
else:
    if 'PYROSETTA_AVAILABLE' in globals() and PYROSETTA_AVAILABLE and pr is not None:
        try:
            pr.init(f'-ignore_unrecognized_res -ignore_zero_occupancy -mute all -holes:dalphaball {advanced_settings["dalphaball_path"]} -corrections::beta_nov16 true -relax:default_repeats 1')
            print("PyRosetta initialized successfully.")
            use_pyrosetta = True
        except Exception as e:
            print(f"PyRosetta detected but failed to initialize: {e}")
            print("Falling back to OpenMM and Biopython routines.")
    else:
        print("PyRosetta not found. Using OpenMM and Biopython routines.")

# Ensure binaries are executable with correct PyRosetta mode
ensure_binaries_executable(use_pyrosetta=use_pyrosetta)

print(f"Running binder design for target {settings_file}")
print(f"Design settings used: {advanced_file}")
print(f"Filtering designs based on {filters_file}")

####################################
# initialise counters
script_start_time = time.time()
trajectory_n = 1
accepted_designs = 0

### start design loop
while True:
    ### check if we have the target number of binders
    final_designs_reached = check_accepted_designs(design_paths, mpnn_csv, final_labels, final_csv, advanced_settings, target_settings, design_labels)

    if final_designs_reached:
        # stop design loop execution
        break

    ### check if we reached maximum allowed trajectories
    max_trajectories_reached = check_n_trajectories(design_paths, advanced_settings)

    if max_trajectories_reached:
        break

    ### Initialise design
    # measure time to generate design
    trajectory_start_time = time.time()

    # generate random seed to vary designs
    seed = int(np.random.randint(0, high=999999, size=1, dtype=int)[0])

    # sample binder design length randomly from defined distribution
    samples = np.arange(min(target_settings["lengths"]), max(target_settings["lengths"]) + 1)
    length = np.random.choice(samples)

    # load desired helicity value to sample different secondary structure contents
    helicity_value = load_helicity(advanced_settings)

    # generate design name and check if same trajectory was already run
    design_name = target_settings["binder_name"] + "_l" + str(length) + "_s"+ str(seed)
    trajectory_dirs = ["Trajectory", "Trajectory/Relaxed", "Trajectory/LowConfidence", "Trajectory/Clashing"]
    trajectory_exists = any(os.path.exists(os.path.join(design_paths[trajectory_dir], design_name + ".pdb")) for trajectory_dir in trajectory_dirs)

    if not trajectory_exists:
        print("Starting trajectory: "+design_name)

        ### Begin binder hallucination
        trajectory = binder_hallucination(design_name, target_settings["starting_pdb"], target_settings["chains"],
                                            target_settings["target_hotspot_residues"], length, seed, helicity_value,
                                            design_models, advanced_settings, design_paths, failure_csv)
        trajectory_metrics = copy_dict(trajectory._tmp["best"]["aux"]["log"]) # contains plddt, ptm, i_ptm, pae, i_pae
        trajectory_pdb = os.path.join(design_paths["Trajectory"], design_name + ".pdb")

        # round the metrics to two decimal places
        trajectory_metrics = {k: round(v, 2) if isinstance(v, float) else v for k, v in trajectory_metrics.items()}

        # time trajectory
        trajectory_time = time.time() - trajectory_start_time
        trajectory_time_text = f"{'%d hours, %d minutes, %d seconds' % (int(trajectory_time // 3600), int((trajectory_time % 3600) // 60), int(trajectory_time % 60))}"
        print("Starting trajectory took: "+trajectory_time_text)
        print("")

        # Proceed if there is no trajectory termination signal
        if trajectory.aux["log"]["terminate"] == "":
            # Relax binder to calculate statistics
            trajectory_relaxed = os.path.join(design_paths["Trajectory/Relaxed"], design_name + ".pdb")
            pr_relax(trajectory_pdb, trajectory_relaxed, use_pyrosetta=use_pyrosetta)

            # define binder chain, placeholder in case multi-chain parsing in ColabDesign gets changed
            binder_chain = "B"

            # Calculate clashes before and after relaxation
            num_clashes_trajectory = calculate_clash_score(trajectory_pdb)
            num_clashes_relaxed = calculate_clash_score(trajectory_relaxed)

            # secondary structure content of starting trajectory binder and interface
            trajectory_alpha, trajectory_beta, trajectory_loops, trajectory_alpha_interface, trajectory_beta_interface, trajectory_loops_interface, trajectory_i_plddt, trajectory_ss_plddt = calc_ss_percentage(trajectory_pdb, advanced_settings, binder_chain)

            # analyze interface scores for relaxed af2 trajectory
            trajectory_interface_scores, trajectory_interface_AA, trajectory_interface_residues = score_interface(trajectory_relaxed, binder_chain, use_pyrosetta=use_pyrosetta)

            # starting binder sequence
            trajectory_sequence = trajectory.get_seq(get_best=True)[0]

            # analyze sequence
            traj_seq_notes = validate_design_sequence(trajectory_sequence, num_clashes_relaxed, advanced_settings)

            # target structure RMSD compared to input PDB
            trajectory_target_rmsd = target_pdb_rmsd(trajectory_pdb, target_settings["starting_pdb"], target_settings["chains"])

            # save trajectory statistics into CSV
            trajectory_data = [design_name, advanced_settings["design_algorithm"], length, seed, helicity_value, target_settings["target_hotspot_residues"], trajectory_sequence, trajectory_interface_residues, 
                                trajectory_metrics['plddt'], trajectory_metrics['ptm'], trajectory_metrics['i_ptm'], trajectory_metrics['pae'], trajectory_metrics['i_pae'],
                                trajectory_i_plddt, trajectory_ss_plddt, num_clashes_trajectory, num_clashes_relaxed, trajectory_interface_scores['binder_score'],
                                trajectory_interface_scores['surface_hydrophobicity'], trajectory_interface_scores['interface_sc'], trajectory_interface_scores['interface_packstat'],
                                trajectory_interface_scores['interface_dG'], trajectory_interface_scores['interface_dSASA'], trajectory_interface_scores['interface_dG_SASA_ratio'],
                                trajectory_interface_scores['interface_fraction'], trajectory_interface_scores['interface_hydrophobicity'], trajectory_interface_scores['interface_nres'], trajectory_interface_scores['interface_interface_hbonds'],
                                trajectory_interface_scores['interface_hbond_percentage'], trajectory_interface_scores['interface_delta_unsat_hbonds'], trajectory_interface_scores['interface_delta_unsat_hbonds_percentage'],
                                trajectory_alpha_interface, trajectory_beta_interface, trajectory_loops_interface, trajectory_alpha, trajectory_beta, trajectory_loops, trajectory_interface_AA, trajectory_target_rmsd, 
                                trajectory_time_text, traj_seq_notes, settings_file, filters_file, advanced_file]
            insert_data(trajectory_csv, trajectory_data)

            # Skip MPNN optimization if no interface residues (no hotspot contact)
            if not trajectory_interface_residues:
                print(f"No interface residues found for {design_name}, skipping MPNN optimization")
                continue
            
            if advanced_settings["enable_mpnn"]:
                # initialise MPNN counters
                mpnn_n = 1
                accepted_mpnn = 0
                mpnn_dict = {}
                design_start_time = time.time()

                ### MPNN redesign of starting binder
                mpnn_trajectories = mpnn_gen_sequence(trajectory_pdb, binder_chain, trajectory_interface_residues, advanced_settings)
                
                existing_mpnn_sequences = set()
                if os.path.exists(mpnn_csv) and os.path.getsize(mpnn_csv) > 0:
                    try:
                        df_mpnn = pd.read_csv(mpnn_csv, usecols=['Sequence'])
                        if not df_mpnn.empty:
                            existing_mpnn_sequences = set(df_mpnn['Sequence'].dropna().astype(str).values)
                    except pd.errors.EmptyDataError:
                        print(f"Warning: {mpnn_csv} is empty or has no columns. Starting with no existing MPNN sequences.")
                    except KeyError:
                        print(f"Warning: 'Sequence' column not found in {mpnn_csv}. Starting with no existing MPNN sequences.")
                    except Exception as e:
                        print(f"Warning: Could not read existing MPNN sequences from {mpnn_csv} due to: {e}. Starting with no existing MPNN sequences.")
                else:
                    print(f"Info: {mpnn_csv} does not exist or is empty. Starting with no existing MPNN sequences.")

                # create set of MPNN sequences with allowed amino acid composition
                restricted_AAs = set(aa.strip().upper() for aa in advanced_settings["omit_AAs"].split(',')) if advanced_settings["force_reject_AA"] else set()

                mpnn_sequences = sorted({
                    mpnn_trajectories['seq'][n][-length:]: {
                        'seq': mpnn_trajectories['seq'][n][-length:],
                        'score': mpnn_trajectories['score'][n],
                        'seqid': mpnn_trajectories['seqid'][n]
                    } for n in range(advanced_settings["num_seqs"])
                    if (not restricted_AAs or not any(aa in mpnn_trajectories['seq'][n][-length:].upper() for aa in restricted_AAs))
                    and mpnn_trajectories['seq'][n][-length:] not in existing_mpnn_sequences
                }.values(), key=lambda x: x['score'])

                del existing_mpnn_sequences
  
                # check whether any sequences are left after amino acid rejection and duplication check, and if yes proceed with prediction
                if mpnn_sequences:
                    # add optimisation for increasing recycles if trajectory is beta sheeted
                    if advanced_settings["optimise_beta"] and float(trajectory_beta) > 15:
                        advanced_settings["num_recycles_validation"] = advanced_settings["optimise_beta_recycles_valid"]

                    ### Compile prediction models once for faster prediction of MPNN sequences
                    clear_mem()
                    # compile complex prediction model
                    complex_prediction_model = mk_afdesign_model(protocol="binder", num_recycles=advanced_settings["num_recycles_validation"], data_dir=advanced_settings["af_params_dir"], 
                                                                use_multimer=multimer_validation, use_initial_guess=advanced_settings["predict_initial_guess"], use_initial_atom_pos=advanced_settings["predict_bigbang"])
                    if advanced_settings["predict_initial_guess"] or advanced_settings["predict_bigbang"]:
                        complex_prediction_model.prep_inputs(pdb_filename=trajectory_pdb, chain='A', binder_chain='B', binder_len=length, use_binder_template=True, rm_target_seq=advanced_settings["rm_template_seq_predict"],
                                                            rm_target_sc=advanced_settings["rm_template_sc_predict"], rm_template_ic=True)
                    else:
                        complex_prediction_model.prep_inputs(pdb_filename=target_settings["starting_pdb"], chain=target_settings["chains"], binder_len=length, rm_target_seq=advanced_settings["rm_template_seq_predict"],
                                                            rm_target_sc=advanced_settings["rm_template_sc_predict"])

                    # compile binder monomer prediction model
                    binder_prediction_model = mk_afdesign_model(protocol="hallucination", use_templates=False, initial_guess=False, 
                                                                use_initial_atom_pos=False, num_recycles=advanced_settings["num_recycles_validation"], 
                                                                data_dir=advanced_settings["af_params_dir"], use_multimer=multimer_validation)
                    binder_prediction_model.prep_inputs(length=length)

                    # iterate over designed sequences        
                    for mpnn_sequence in mpnn_sequences:
                        mpnn_time = time.time()

                        # generate mpnn design name numbering
                        mpnn_design_name = design_name + "_mpnn" + str(mpnn_n)
                        mpnn_score = round(mpnn_sequence['score'],2)
                        mpnn_seqid = round(mpnn_sequence['seqid'],2)

                        # add design to dictionary
                        mpnn_dict[mpnn_design_name] = {'seq': mpnn_sequence['seq'], 'score': mpnn_score, 'seqid': mpnn_seqid}

                        # save fasta sequence
                        if advanced_settings["save_mpnn_fasta"] is True:
                            save_fasta(mpnn_design_name, mpnn_sequence['seq'], design_paths)
                        
                        ### Predict mpnn redesigned binder complex using masked templates
                        mpnn_complex_statistics, pass_af2_filters, early_filter_details = predict_binder_complex(complex_prediction_model,
                                                                                        mpnn_sequence['seq'], mpnn_design_name,
                                                                                        target_settings["starting_pdb"], target_settings["chains"],
                                                                                        length, trajectory_pdb, prediction_models, advanced_settings,
                                                                                        filters, design_paths, failure_csv, use_pyrosetta=use_pyrosetta)

                        # if AF2 filters are not passed then skip the scoring but log the failure
                        if not pass_af2_filters:
                            print(f"Base AF2 filters not passed for {mpnn_design_name}, skipping full interface scoring.")

                            # Log to rejected_mpnn_full_stats.csv for early AF2 failures
                            rejected_data_list_for_csv = [mpnn_design_name, mpnn_sequence['seq']]
                            
                            failed_base_metrics_for_this_design = set()
                            if early_filter_details: # early_filter_details is the filter_failures dict e.g. {"1_pLDDT": 1}
                                for specific_model_failure_key in early_filter_details.keys():
                                    parts = specific_model_failure_key.split('_')
                                    if len(parts) > 1 and parts[0].isdigit(): # e.g. "1_pLDDT" -> "pLDDT"
                                        base_metric_name = ''.join(parts[1:]) # Corrected: was '_'.join, should be '' to match pLDDT from pLDDT
                                        # Special case for i_pTM, i_pAE, i_pLDDT as they already contain an underscore
                                        if parts[1] == "i" and len(parts) > 2: # e.g. "1_i_pTM" -> "i_pTM"
                                            base_metric_name = parts[1] + "_" + ''.join(parts[2:])
                                        failed_base_metrics_for_this_design.add(base_metric_name)
                                    else: # Should not happen with current predict_binder_complex failures, but good for robustness
                                        failed_base_metrics_for_this_design.add(specific_model_failure_key)

                            for base_filter_col_name_in_log_csv in filter_column_names_for_rejected_log:
                                if base_filter_col_name_in_log_csv in failed_base_metrics_for_this_design:
                                    rejected_data_list_for_csv.append(1)
                                else:
                                    rejected_data_list_for_csv.append(0)
                            insert_data(rejected_mpnn_full_stats_csv, rejected_data_list_for_csv)
                            
                            mpnn_n += 1
                            continue

                        # calculate statistics for each model individually
                        for model_num in prediction_models:
                            mpnn_design_pdb = os.path.join(design_paths["MPNN"], f"{mpnn_design_name}_model{model_num+1}.pdb")
                            mpnn_design_relaxed = os.path.join(design_paths["MPNN/Relaxed"], f"{mpnn_design_name}_model{model_num+1}.pdb")

                            if os.path.exists(mpnn_design_pdb):
                                # Calculate clashes before and after relaxation
                                num_clashes_mpnn = calculate_clash_score(mpnn_design_pdb)
                                num_clashes_mpnn_relaxed = calculate_clash_score(mpnn_design_relaxed)

                                # analyze interface scores for relaxed af2 trajectory
                                mpnn_interface_scores, mpnn_interface_AA, mpnn_interface_residues = score_interface(mpnn_design_relaxed, binder_chain, use_pyrosetta=use_pyrosetta)

                                # secondary structure content of starting trajectory binder
                                mpnn_alpha, mpnn_beta, mpnn_loops, mpnn_alpha_interface, mpnn_beta_interface, mpnn_loops_interface, mpnn_i_plddt, mpnn_ss_plddt = calc_ss_percentage(mpnn_design_pdb, advanced_settings, binder_chain)
                                
                                # unaligned RMSD calculate to determine if binder is in the designed binding site
                                rmsd_site = unaligned_rmsd(trajectory_pdb, mpnn_design_pdb, binder_chain, binder_chain, use_pyrosetta=use_pyrosetta)

                                # calculate RMSD of target compared to input PDB
                                target_rmsd = target_pdb_rmsd(mpnn_design_pdb, target_settings["starting_pdb"], target_settings["chains"])

                                # add the additional statistics to the mpnn_complex_statistics dictionary
                                mpnn_complex_statistics[model_num+1].update({
                                    'i_pLDDT': mpnn_i_plddt,
                                    'ss_pLDDT': mpnn_ss_plddt,
                                    'Unrelaxed_Clashes': num_clashes_mpnn,
                                    'Relaxed_Clashes': num_clashes_mpnn_relaxed,
                                    'Binder_Energy_Score': mpnn_interface_scores['binder_score'],
                                    'Surface_Hydrophobicity': mpnn_interface_scores['surface_hydrophobicity'],
                                    'ShapeComplementarity': mpnn_interface_scores['interface_sc'],
                                    'PackStat': mpnn_interface_scores['interface_packstat'],
                                    'dG': mpnn_interface_scores['interface_dG'],
                                    'dSASA': mpnn_interface_scores['interface_dSASA'], 
                                    'dG/dSASA': mpnn_interface_scores['interface_dG_SASA_ratio'],
                                    'Interface_SASA_%': mpnn_interface_scores['interface_fraction'],
                                    'Interface_Hydrophobicity': mpnn_interface_scores['interface_hydrophobicity'],
                                    'n_InterfaceResidues': mpnn_interface_scores['interface_nres'],
                                    'n_InterfaceHbonds': mpnn_interface_scores['interface_interface_hbonds'],
                                    'InterfaceHbondsPercentage': mpnn_interface_scores['interface_hbond_percentage'],
                                    'n_InterfaceUnsatHbonds': mpnn_interface_scores['interface_delta_unsat_hbonds'],
                                    'InterfaceUnsatHbondsPercentage': mpnn_interface_scores['interface_delta_unsat_hbonds_percentage'],
                                    'InterfaceAAs': mpnn_interface_AA,
                                    'Interface_Helix%': mpnn_alpha_interface,
                                    'Interface_BetaSheet%': mpnn_beta_interface,
                                    'Interface_Loop%': mpnn_loops_interface,
                                    'Binder_Helix%': mpnn_alpha,
                                    'Binder_BetaSheet%': mpnn_beta,
                                    'Binder_Loop%': mpnn_loops,
                                    'Hotspot_RMSD': rmsd_site,
                                    'Target_RMSD': target_rmsd
                                })

                                # save space by removing unrelaxed predicted mpnn complex pdb?
                                if advanced_settings["remove_unrelaxed_complex"]:
                                    os.remove(mpnn_design_pdb)

                        # calculate complex averages
                        mpnn_complex_averages = calculate_averages(mpnn_complex_statistics, handle_aa=True)
                        
                        ### Predict binder alone in single sequence mode
                        binder_statistics = predict_binder_alone(binder_prediction_model, mpnn_sequence['seq'], mpnn_design_name, length,
                                                                trajectory_pdb, binder_chain, prediction_models, advanced_settings, design_paths, 
                                                                use_pyrosetta=use_pyrosetta)

                        # extract RMSDs of binder to the original trajectory
                        for model_num in prediction_models:
                            mpnn_binder_pdb = os.path.join(design_paths["MPNN/Binder"], f"{mpnn_design_name}_model{model_num+1}.pdb")

                            if os.path.exists(mpnn_binder_pdb):
                                rmsd_binder = unaligned_rmsd(trajectory_pdb, mpnn_binder_pdb, binder_chain, "A", use_pyrosetta=use_pyrosetta)

                            # append to statistics
                            binder_statistics[model_num+1].update({
                                    'Binder_RMSD': rmsd_binder
                                })

                            # save space by removing binder monomer models?
                            if advanced_settings["remove_binder_monomer"]:
                                os.remove(mpnn_binder_pdb)

                        # calculate binder averages
                        binder_averages = calculate_averages(binder_statistics)

                        # analyze sequence to make sure there are no cysteins and it contains residues that absorb UV for detection
                        seq_notes = validate_design_sequence(mpnn_sequence['seq'], mpnn_complex_averages.get('Relaxed_Clashes', None), advanced_settings)

                        # measure time to generate design
                        mpnn_end_time = time.time() - mpnn_time
                        elapsed_mpnn_text = f"{'%d hours, %d minutes, %d seconds' % (int(mpnn_end_time // 3600), int((mpnn_end_time % 3600) // 60), int(mpnn_end_time % 60))}"


                        # Insert statistics about MPNN design into CSV, will return None if corresponding model does note exist
                        model_numbers = range(1, 6)
                        statistics_labels = ['pLDDT', 'pTM', 'i_pTM', 'pAE', 'i_pAE', 'i_pLDDT', 'ss_pLDDT', 'Unrelaxed_Clashes', 'Relaxed_Clashes', 'Binder_Energy_Score', 'Surface_Hydrophobicity',
                                            'ShapeComplementarity', 'PackStat', 'dG', 'dSASA', 'dG/dSASA', 'Interface_SASA_%', 'Interface_Hydrophobicity', 'n_InterfaceResidues', 'n_InterfaceHbonds', 'InterfaceHbondsPercentage',
                                            'n_InterfaceUnsatHbonds', 'InterfaceUnsatHbondsPercentage', 'Interface_Helix%', 'Interface_BetaSheet%', 'Interface_Loop%', 'Binder_Helix%',
                                            'Binder_BetaSheet%', 'Binder_Loop%', 'InterfaceAAs', 'Hotspot_RMSD', 'Target_RMSD']

                        # Initialize mpnn_data with the non-statistical data
                        mpnn_data = [mpnn_design_name, advanced_settings["design_algorithm"], length, seed, helicity_value, target_settings["target_hotspot_residues"], mpnn_sequence['seq'], mpnn_interface_residues, mpnn_score, mpnn_seqid]

                        # Add the statistical data for mpnn_complex
                        for label in statistics_labels:
                            mpnn_data.append(mpnn_complex_averages.get(label, None))
                            for model in model_numbers:
                                mpnn_data.append(mpnn_complex_statistics.get(model, {}).get(label, None))

                        # Add the statistical data for binder
                        for label in ['pLDDT', 'pTM', 'pAE', 'Binder_RMSD']:  # These are the labels for binder alone
                            mpnn_data.append(binder_averages.get(label, None))
                            for model in model_numbers:
                                mpnn_data.append(binder_statistics.get(model, {}).get(label, None))

                        # Add the remaining non-statistical data
                        mpnn_data.extend([elapsed_mpnn_text, seq_notes, settings_file, filters_file, advanced_file])

                        # insert data into csv
                        insert_data(mpnn_csv, mpnn_data)

                        # find best model number by pLDDT
                        plddt_values = {i: mpnn_data[i] for i in range(11, 15) if mpnn_data[i] is not None}

                        # Find the key with the highest value
                        highest_plddt_key = int(max(plddt_values, key=plddt_values.get))

                        # Output the number part of the key
                        best_model_number = highest_plddt_key - 10
                        best_model_pdb = os.path.join(design_paths["MPNN/Relaxed"], f"{mpnn_design_name}_model{best_model_number}.pdb")

                        # run design data against filter thresholds
                        filter_conditions = check_filters(mpnn_data, design_labels, filters)
                        if filter_conditions == True:
                            print(mpnn_design_name+" passed all filters")
                            accepted_mpnn += 1
                            accepted_designs += 1
                            
                            # copy designs to accepted folder
                            shutil.copy(best_model_pdb, design_paths["Accepted"])

                            # insert data into final csv
                            final_data = [''] + mpnn_data
                            insert_data(final_csv, final_data)

                            # copy animation from accepted trajectory
                            if advanced_settings["save_design_animations"]:
                                accepted_animation = os.path.join(design_paths["Accepted/Animation"], f"{design_name}.html")
                                if not os.path.exists(accepted_animation):
                                    shutil.copy(os.path.join(design_paths["Trajectory/Animation"], f"{design_name}.html"), accepted_animation)

                            # copy plots of accepted trajectory
                            plot_files = os.listdir(design_paths["Trajectory/Plots"])
                            plots_to_copy = [f for f in plot_files if f.startswith(design_name) and f.endswith('.png')]
                            for accepted_plot in plots_to_copy:
                                source_plot = os.path.join(design_paths["Trajectory/Plots"], accepted_plot)
                                target_plot = os.path.join(design_paths["Accepted/Plots"], accepted_plot)
                                if not os.path.exists(target_plot):
                                    shutil.copy(source_plot, target_plot)

                        else:
                            print(f"Unmet filter conditions for {mpnn_design_name}")
                            failure_df = pd.read_csv(failure_csv)
                            special_prefixes = ('Average_', '1_', '2_', '3_', '4_', '5_')
                            incremented_columns = set()

                            for column in filter_conditions:
                                base_column = column
                                for prefix in special_prefixes:
                                    if column.startswith(prefix):
                                        base_column = column.split('_', 1)[1]
                                        break # Corrected: was missing break

                                if base_column not in incremented_columns:
                                    if base_column in failure_df.columns: # Check if column exists before incrementing
                                        failure_df[base_column] = failure_df[base_column] + 1
                                    else:
                                        # This case should ideally not happen if generate_filter_pass_csv creates all potential base_columns
                                        print(f"Warning: Base column '{base_column}' not found in {failure_csv}. It won't be incremented for this failure.")
                                    incremented_columns.add(base_column)
                            
                            failure_df.to_csv(failure_csv, index=False)

                            # Log to rejected_mpnn_full_stats.csv
                            rejected_data_list_for_csv = [mpnn_design_name, mpnn_sequence['seq']]
                            for filter_col_name in filter_column_names_for_rejected_log:
                                if filter_col_name in incremented_columns:
                                    rejected_data_list_for_csv.append(1)
                                else:
                                    rejected_data_list_for_csv.append(0)
                            insert_data(rejected_mpnn_full_stats_csv, rejected_data_list_for_csv)
                            
                            shutil.copy(best_model_pdb, design_paths["Rejected"])
                        
                        # increase MPNN design number
                        mpnn_n += 1
                        
                        # Force garbage collection after each MPNN design to prevent file descriptor accumulation
                        gc.collect()

                        # if enough mpnn sequences of the same trajectory pass filters then stop
                        if accepted_mpnn >= advanced_settings["max_mpnn_sequences"]:
                            break

                    if accepted_mpnn >= 1:
                        print("Found "+str(accepted_mpnn)+" MPNN designs passing filters")
                        print("")
                    else:
                        print("No accepted MPNN designs found for this trajectory.")
                        print("")

                else:
                    print('Duplicate MPNN designs sampled with different trajectory, skipping current trajectory optimisation')
                    print("")

                # save space by removing unrelaxed design trajectory PDB
                if advanced_settings["remove_unrelaxed_trajectory"]:
                    os.remove(trajectory_pdb)

                # measure time it took to generate designs for one trajectory
                design_time = time.time() - design_start_time
                design_time_text = f"{'%d hours, %d minutes, %d seconds' % (int(design_time // 3600), int((design_time % 3600) // 60), int(design_time % 60))}"
                print("Design and validation of trajectory "+design_name+" took: "+design_time_text)

            # analyse the rejection rate of trajectories to see if we need to readjust the design weights
            if trajectory_n >= advanced_settings["start_monitoring"] and advanced_settings["enable_rejection_check"]:
                acceptance = accepted_designs / trajectory_n
                if not acceptance >= advanced_settings["acceptance_rate"]:
                    print("The ratio of successful designs is lower than defined acceptance rate! Consider changing your design settings!")
                    print("Script execution stopping...")
                    break

        # increase trajectory number
        trajectory_n += 1
        
        # Force garbage collection and clear DSSP cache every 10 trajectories to prevent memory/fd accumulation
        if trajectory_n % 10 == 0:
            clear_dssp_cache()
            print(f"Cleared DSSP cache after {trajectory_n} trajectories")
        
        # Force garbage collection more frequently to prevent file descriptor accumulation
        gc.collect()

### Script finished
elapsed_time = time.time() - script_start_time
elapsed_text = f"{'%d hours, %d minutes, %d seconds' % (int(elapsed_time // 3600), int((elapsed_time % 3600) // 60), int(elapsed_time % 60))}"
print("Finished all designs. Script execution for "+str(trajectory_n)+" trajectories took: "+elapsed_text)