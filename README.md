# BindCraft (v1.5) with Optional PyRosetta Bypass

This repository contains a modified version of Martin Pacesa's BindCraft (v1.5). The primary change is the introduction of an **optional PyRosetta bypass mechanism**.

For comprehensive details on the original BindCraft pipeline, features, advanced settings, and filter explanations, please refer to the **original BindCraft repository: [https://github.com/martinpacesa/BindCraft](https://github.com/martinpacesa/BindCraft)** and the [original preprint](https://www.biorxiv.org/content/10.1101/2024.09.30.615802).

## Key Modification: PyRosetta Bypass Functionality

The `--no-pyrosetta` flag, usable during both installation and runtime, alters BindCraft's behavior:

*   **OpenMM Relaxation:** When PyRosetta is bypassed, structural relaxation is performed using an OpenMM-based protocol instead of PyRosetta's FastRelax. This process includes structure preparation with PDBFixer, ramped backbone restraints to prevent excessive backbone movement, implicit solvation (OBC2), a custom repulsive force to mitigate clashes, and MD shakes.
*   **Interface Scoring Bypass:** PyRosetta-dependent interface scoring (e.g., dG, ShapeComplementarity) is skipped. Placeholder values are returned for these metrics, configured to pass default filter thresholds.
*   **Biopython Alternatives:** Critical functions normally reliant on PyRosetta (e.g., for RMSD calculations, PDB alignment) are substituted with Biopython-based implementations.

This bypass is intended for situations where:
*   PyRosetta licensing or installation poses a challenge.
*   Speed is prioritized over detailed Rosetta-based structural metrics.
*   A reduced dependency footprint is desired.

**Important Note:** When using the PyRosetta bypass, structural analysis and filtering based on Rosetta-specific metrics are not performed. Evaluate design quality accordingly.

## Installation

1.  Clone this modified repository:
    ```bash
    git clone https://github.com/cytokineking/PyRosetta-Optional-BindCraft-v1.5 [install_folder]
    ```
2.  Navigate into your install folder (`cd [install_folder]`) and run the installation script. A CUDA-compatible Nvidia graphics card is required.

    Installation script options:
    *   `--cuda CUDAVERSION`: Specify your CUDA version (e.g., '12.4').
    *   `--pkg_manager MANAGER`: Specify 'mamba' or 'conda' (default: 'conda').
    *   `--no-pyrosetta`: **Use this flag to install without PyRosetta dependencies and enable the bypass mode features.**

    **Example: PyRosetta-Free Installation (Recommended for bypass)**
    ```bash
    bash install_bindcraft.sh --cuda '12.4' --pkg_manager 'conda' --no-pyrosetta
    ```

    **Example: Standard Installation (Includes PyRosetta)**
    If you choose to install with PyRosetta, a license is required for commercial use.
    ```bash
    bash install_bindcraft.sh --cuda '12.4' --pkg_manager 'conda'
    ```

## Running BindCraft

Activate the Conda environment:
```bash
conda activate BindCraft
```

Change to your BindCraft installation directory:
```bash
# cd /path/to/your/bindcraft/folder/
```

**To use the PyRosetta bypass at runtime, include the `--no-pyrosetta` flag:**
```bash
python -u ./bindcraft.py --settings './settings_target/your_target.json' --filters './settings_filters/default_filters.json' --advanced './settings_advanced/default_4stage_multimer.json' --no-pyrosetta
```

**To run with PyRosetta enabled (assuming it was installed):**
```bash
python -u ./bindcraft.py --settings './settings_target/your_target.json' --filters './settings_filters/default_filters.json' --advanced './settings_advanced/default_4stage_multimer.json'
```

**Note:** Even if you installed BindCraft *with* PyRosetta, you can still run in bypass mode by adding the `--no-pyrosetta` flag at runtime.

For details on configuring target settings (`--settings`), filters (`--filters`), advanced parameters (`--advanced`), and other operational aspects of BindCraft, please consult the documentation in the [original BindCraft repository](https://github.com/martinpacesa/BindCraft).
