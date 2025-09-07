import argparse
import os
import sys
import json
from datetime import datetime

# Ensure the project root is in the Python path to allow imports from 'functions'
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

try:
    from functions.pyrosetta_utils import pr_relax, PYROSETTA_AVAILABLE
    from functions.pr_alternative_utils import openmm_relax
    from functions.generic_utils import clean_pdb
    if PYROSETTA_AVAILABLE: # Import pr here if available, for initialization
        from pyrosetta import init as pr_init
except ImportError as e:
    print(f"Critical Import Error: {e}")
    print(f"Attempted to add project root '{project_root}' to sys.path.")
    print(f"Current sys.path: {sys.path}")
    print("Please ensure that the 'functions' directory is a Python package (contains __init__.py)")
    print("and that all dependencies of 'functions.pyrosetta_utils.openmm_relax' are installed (e.g., openmm, pdbfixer).")
    sys.exit(1)
except Exception as e_gen:
    print(f"An unexpected error occurred during imports: {e_gen}")
    sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Test OpenMM and PyRosetta relaxation with optional performance sweeps. "
            "Prints per-stage timing and energy metrics."
        )
    )
    parser.add_argument(
        "input_pdb_path",
        type=str,
        help="Path to the input PDB file for relaxation."
    )
    parser.add_argument(
        "base_output_pdb_path",
        type=str,
        help="Base path for the relaxed output PDB files. Suffixes _openmm.pdb and _pyrosetta.pdb will be added."
    )

    # OpenMM performance controls (single-run)
    parser.add_argument("--platform", choices=["auto", "cpu", "gpu", "cuda", "opencl"], default="auto",
                        help="Platform selection for OpenMM relax (auto uses default search)")
    parser.add_argument("--md-steps-per-shake", type=int, default=5000,
                        help="MD steps per shake (0 disables MD shakes)")
    parser.add_argument("--restraint-ramp", type=str, default="1.0,0.4,0.0",
                        help="Comma-separated restraint ramp factors (e.g., 1.0,0.4,0.0)")
    parser.add_argument("--lj-ramp", type=str, default="0.0,1.5,3.0",
                        help="Comma-separated LJ repulsion ramp factors (e.g., 0.0,1.5,3.0)")
    parser.add_argument("--max-iters", type=int, default=1000,
                        help="Max iterations per minimization stage (0 for unlimited)")
    parser.add_argument("--ramp-force-tol", type=float, default=2.0,
                        help="Force tolerance (kJ/mol/nm) for ramp minimization stages")
    parser.add_argument("--final-force-tol", type=float, default=0.1,
                        help="Force tolerance (kJ/mol/nm) for final stage minimization")
    parser.add_argument("--json-out", type=str, default=None,
                        help="Optional path to write single-run perf JSON")

    # FASPR repacking toggle
    parser.add_argument("--faspr", action="store_true",
                        help="Run FASPR side-chain repacking after OpenMM relax (requires functions/FASPR and rotamer lib)")

    # Sweep controls
    parser.add_argument("--sweep-basic", action="store_true",
                        help="Run a basic grid sweep over CPU/GPU x MD on/off x ramp stages 1/2/3")
    parser.add_argument("--json-dir", type=str, default=None,
                        help="If set with --sweep-basic, writes one JSON per configuration to this directory")

    args = parser.parse_args()

    print(f"--- Test Script: Starting Relaxations ---")
    print(f"Input PDB: {args.input_pdb_path}")
    print(f"Base Output PDB Path: {args.base_output_pdb_path}")

    if not os.path.exists(args.input_pdb_path):
        print(f"Error: Input PDB file not found: {args.input_pdb_path}")
        sys.exit(1)

    # Clean the input PDB before attempting relaxation
    print(f"Cleaning input PDB file: {args.input_pdb_path}...")
    try:
        clean_pdb(args.input_pdb_path)
        print(f"Input PDB file cleaned successfully.")
    except Exception as e_clean:
        print(f"Error during clean_pdb for {args.input_pdb_path}: {e_clean}")
        print("Proceeding with relaxation using the original (uncleaned) PDB.")

    # Derive specific output paths
    output_dir = os.path.dirname(args.base_output_pdb_path)
    base_name = os.path.basename(args.base_output_pdb_path)
    # root will be the base_name without its original extension (if any)
    root = os.path.splitext(base_name)[0]

    openmm_output_pdb_path = os.path.join(output_dir, root + "_openmm.pdb")
    pyrosetta_output_pdb_path = os.path.join(output_dir, root + "_pyrosetta.pdb")
    
    print(f"Target OpenMM output: {openmm_output_pdb_path}")
    print(f"Target PyRosetta output: {pyrosetta_output_pdb_path}")

    # Ensure the output directory exists
    # The directory for base_output_pdb_path is the same for suffixed files
    if output_dir: 
        try:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
                print(f"Created output directory: {output_dir}")
        except OSError as e:
            print(f"Error creating output directory {output_dir}: {e}")
            sys.exit(1)
            
    # --- OpenMM Relaxation ---
    def _parse_floats_csv(s, default):
        try:
            vals = [float(x.strip()) for x in str(s).split(',') if x.strip() != ""]
            return vals if vals else list(default)
        except Exception:
            return list(default)

    def _platform_override_from_choice(choice):
        c = str(choice).lower()
        if c == "cpu":
            return ["CPU"], False
        if c == "cuda":
            return ["CUDA"], True
        if c == "opencl":
            return ["OpenCL"], True
        if c == "gpu":
            return ["OpenCL", "CUDA"], True
        return None, True  # auto

    def _run_openmm_once(tag, input_pdb, out_pdb, platform_choice, md_steps, rest_ramp, lj_ramp, max_iters, ramp_tol, final_tol, json_out_path=None, faspr=False):
        print(f"\n--- OpenMM Relax Run: {tag} ---")
        perf = {}
        override, use_gpu = _platform_override_from_choice(platform_choice)
        try:
            used_platform = openmm_relax(
                input_pdb,
                out_pdb,
                use_gpu_relax=use_gpu,
                openmm_max_iterations=max_iters,
                openmm_ramp_force_tolerance_kj_mol_nm=ramp_tol,
                openmm_final_force_tolerance_kj_mol_nm=final_tol,
                restraint_k_kcal_mol_A2=3.0,
                restraint_ramp_factors=rest_ramp,
                md_steps_per_shake=md_steps,
                lj_rep_base_k_kj_mol=10.0,
                lj_rep_ramp_factors=lj_ramp,
                perf_report=perf,
                override_platform_order=override,
                use_faspr_repack=bool(faspr),
                post_faspr_minimize=True
            )
            print(f"Platform: {used_platform}")
            print(f"Total seconds: {perf.get('total_seconds'):.2f}")
            print(f"Initial minimization seconds: {perf.get('initial_min_seconds', 0.0):.2f}")
            print(f"Ramp count: {perf.get('ramp_count')}, MD steps/shake: {perf.get('md_steps_per_shake')}")
            if perf.get('best_energy_kj') is not None:
                print(f"Best energy (kJ/mol): {perf['best_energy_kj']:.2f}")
            for st in perf.get('stages', []):
                print(
                    f"  Stage {st['stage_index']}: md_steps={st['md_steps_run']}, "
                    f"md_s={st['md_seconds']:.2f}, min_calls={st['min_calls']}, min_s={st['min_seconds']:.2f}, "
                    f"E0={st.get('energy_start_kj')}, Emd={st.get('md_post_energy_kj')}, Efin={st.get('final_energy_kj')}"
                )
            if faspr:
                print(f"FASPR: success={perf.get('faspr_success')}, faspr_s={perf.get('faspr_seconds')}, post_min_s={perf.get('post_faspr_min_seconds')}")
            if json_out_path:
                try:
                    os.makedirs(os.path.dirname(json_out_path), exist_ok=True) if os.path.dirname(json_out_path) else None
                    with open(json_out_path, 'w') as jf:
                        json.dump(perf, jf, indent=2)
                    print(f"Perf JSON written: {json_out_path}")
                except Exception as je:
                    print(f"Failed to write JSON to {json_out_path}: {je}")
            if os.path.exists(out_pdb):
                print(f"OpenMM Relaxed PDB saved to: {out_pdb}")
        except Exception as e:
            print(f"OpenMM run failed for tag '{tag}': {type(e).__name__}: {e}")
            import traceback
            traceback.print_exc()

    # Single-run or sweep execution
    if args.sweep_basic:
        print("\nRunning basic sweep over platform x MD on/off x ramp stages 1/2/3...")
        base_rest = _parse_floats_csv(args.restraint_ramp, [1.0, 0.4, 0.0])
        base_lj = _parse_floats_csv(args.lj_ramp, [0.0, 1.5, 3.0])
        platforms = ["cpu", "gpu"]
        md_opts = [0, max(0, int(args.md_steps_per_shake))]
        stage_counts = [1, 2, 3]
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        for plat in platforms:
            for md in md_opts:
                for sc in stage_counts:
                    rest = base_rest[:sc]
                    lj = base_lj[:sc]
                    tag = f"{plat}_md{md}_st{sc}"
                    out_pdb = os.path.splitext(openmm_output_pdb_path)[0] + f"_{tag}.pdb"
                    json_out = None
                    if args.json_dir:
                        json_out = os.path.join(args.json_dir, f"perf_{tag}_{stamp}.json")
                    _run_openmm_once(
                        tag, args.input_pdb_path, out_pdb,
                        plat, md, rest, lj,
                        args.max_iters, args.ramp_force_tol, args.final_force_tol,
                        json_out_path=json_out,
                        faspr=args.faspr
                    )
    else:
        # Single-run
        rest = _parse_floats_csv(args.restraint_ramp, [1.0, 0.4, 0.0])
        lj = _parse_floats_csv(args.lj_ramp, [0.0, 1.5, 3.0])
        _run_openmm_once(
            "single",
            args.input_pdb_path,
            openmm_output_pdb_path,
            args.platform,
            max(0, int(args.md_steps_per_shake)),
            rest,
            lj,
            args.max_iters,
            args.ramp_force_tol,
            args.final_force_tol,
            json_out_path=args.json_out,
            faspr=args.faspr
        )

    # --- PyRosetta Relaxation ---
    print(f"\n--- Starting PyRosetta Relaxation ---")
    if PYROSETTA_AVAILABLE:
        print("Initializing PyRosetta...")
        try:
            pr_init("-mute all") # Basic initialization for the test script
            print("PyRosetta initialized successfully for the test script.")
        except Exception as e_init:
            print(f"Error during PyRosetta initialization: {e_init}")
            print("Skipping PyRosetta relaxation due to initialization error.")
        else:
            print(f"PyRosetta is available. Calling pr_relax function for {pyrosetta_output_pdb_path}...")
            try:
                # Call pr_relax with use_pyrosetta=True to attempt PyRosetta relaxation
                pr_relax(args.input_pdb_path, pyrosetta_output_pdb_path, use_pyrosetta=True)
                print(f"pr_relax function completed.")
                if os.path.exists(pyrosetta_output_pdb_path):
                    print(f"PyRosetta Relaxed PDB saved to: {pyrosetta_output_pdb_path}")
                else:
                    print(f"Warning: PyRosetta Output PDB file was not created at {pyrosetta_output_pdb_path} by pr_relax.")
            except Exception as e:
                print(f"--- Test Script: Error during pr_relax (PyRosetta) execution ---")
                print(f"An exception occurred: {type(e).__name__}: {e}")
                import traceback
                traceback.print_exc()
    else:
        print(f"PyRosetta is not available (PYROSETTA_AVAILABLE=False). Skipping PyRosetta relaxation.")
        print(f"If you intended to test PyRosetta relaxation, please ensure it's correctly installed and configured.")

    print(f"\n--- Test Script: Finished ---") 