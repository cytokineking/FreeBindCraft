import argparse
import os
import sys

# This script should be in the root of the BindCraft-main-v1.5 project.
# Ensure the project root is in the Python path to allow imports from 'functions'
project_root = os.path.dirname(os.path.abspath(__file__))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

try:
    from functions.pyrosetta_utils import openmm_relax, pr_relax, PYROSETTA_AVAILABLE
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
        description="Test OpenMM and PyRosetta relaxation. "\
                    "Prints verbose output from the relaxation functions."
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
    print(f"\n--- Starting OpenMM Relaxation ---")
    print(f"Calling openmm_relax function for {openmm_output_pdb_path}...")
    try:
        # Call the imported openmm_relax function
        openmm_relax(args.input_pdb_path, openmm_output_pdb_path, use_gpu_relax=True) # Attempt GPU relaxation (CUDA -> OpenCL -> CPU fallback)
        print(f"openmm_relax function completed.")
        if os.path.exists(openmm_output_pdb_path):
            print(f"OpenMM Relaxed PDB saved to: {openmm_output_pdb_path}")
        else:
            print(f"Warning: OpenMM Output PDB file was not created at {openmm_output_pdb_path}.")

    except Exception as e:
        print(f"--- Test Script: Error during openmm_relax execution (forced CPU) ---")
        print(f"An exception occurred: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()

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