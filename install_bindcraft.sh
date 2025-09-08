#!/bin/bash
################## BindCraft installation script
################## specify conda/mamba folder, and installation folder for git repositories, and whether to use mamba or $pkg_manager
# Default value for pkg_manager
pkg_manager='conda'
cuda=''
install_pyrosetta=true

# Define the short and long options
OPTIONS=p:c:n
LONGOPTIONS=pkg_manager:,cuda:,no-pyrosetta

# Parse the command-line options
PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTIONS --name "$0" -- "$@")
eval set -- "$PARSED"

# Process the command-line options
while true; do
  case "$1" in
    -p|--pkg_manager)
      pkg_manager="$2"
      shift 2
      ;;
    -c|--cuda)
      cuda="$2"
      shift 2
      ;;
    -n|--no-pyrosetta)
      install_pyrosetta=false
      shift
      ;;
    --)
      shift
      break
      ;;
    *)
      echo -e "Invalid option $1" >&2
      exit 1
      ;;
  esac
done

# Example usage of the parsed variables
echo -e "Package manager: $pkg_manager"
echo -e "CUDA: $cuda"
echo -e "Install PyRosetta: $install_pyrosetta"

############################################################################################################
############################################################################################################
################## initialisation
SECONDS=0

# set paths needed for installation and check for conda installation
install_dir=$(pwd)
CONDA_BASE=$(conda info --base 2>/dev/null) || { echo -e "Error: conda is not installed or cannot be initialised."; exit 1; }
echo -e "Conda is installed at: $CONDA_BASE"

### BindCraft install begin, create base environment
echo -e "Installing BindCraft environment\n"
$pkg_manager create --name BindCraft python=3.10 -y || { echo -e "Error: Failed to create BindCraft conda environment"; exit 1; }
conda env list | grep -w 'BindCraft' >/dev/null 2>&1 || { echo -e "Error: Conda environment 'BindCraft' does not exist after creation."; exit 1; }

# Load newly created BindCraft environment
echo -e "Loading BindCraft environment\n"
source ${CONDA_BASE}/bin/activate ${CONDA_BASE}/envs/BindCraft || { echo -e "Error: Failed to activate the BindCraft environment."; exit 1; }
[ "$CONDA_DEFAULT_ENV" = "BindCraft" ] || { echo -e "Error: The BindCraft environment is not active."; exit 1; }
echo -e "BindCraft environment activated at ${CONDA_BASE}/envs/BindCraft"

# install required conda packages
echo -e "Installing conda requirements\n"

# Base packages (needed regardless of PyRosetta)
BASE_PACKAGES="pip pandas matplotlib numpy<2.0.0 biopython scipy pdbfixer openmm seaborn libgfortran5 tqdm jupyter ffmpeg fsspec py3dmol chex dm-haiku flax<0.10.0 dm-tree joblib ml-collections immutabledict optax"

# Pin JAX to stable version 0.6.0
echo -e "Using JAX/jaxlib version 0.6.0 for stability\n"

# Install packages with or without PyRosetta
if [ "$install_pyrosetta" = true ]; then
    echo -e "Installing with PyRosetta\n"
    if [ -n "$cuda" ]; then
        CONDA_OVERRIDE_CUDA="$cuda" $pkg_manager install \
            $BASE_PACKAGES pyrosetta "jaxlib=0.6.0=*cuda*" "jax=0.6.0" cuda-nvcc cudnn \
            -c conda-forge -c nvidia -c "https://conda.rosettacommons.org" -y \
            || { echo -e "Error: Failed to install conda packages with PyRosetta."; exit 1; }
    else
        $pkg_manager install \
            $BASE_PACKAGES pyrosetta "jaxlib=0.6.0" "jax=0.6.0" cuda-nvcc cudnn \
            -c conda-forge -c nvidia -c "https://conda.rosettacommons.org" -y \
            || { echo -e "Error: Failed to install conda packages with PyRosetta."; exit 1; }
    fi
else
    echo -e "Installing without PyRosetta\n"
    if [ -n "$cuda" ]; then
        CONDA_OVERRIDE_CUDA="$cuda" $pkg_manager install \
            $BASE_PACKAGES "jaxlib=0.6.0=*cuda*" "jax=0.6.0" cuda-nvcc cudnn \
            -c conda-forge -c nvidia -y \
            || { echo -e "Error: Failed to install conda packages without PyRosetta."; exit 1; }
    else
        $pkg_manager install \
            $BASE_PACKAGES "jaxlib=0.6.0" "jax=0.6.0" cuda-nvcc cudnn \
            -c conda-forge -c nvidia -y \
            || { echo -e "Error: Failed to install conda packages without PyRosetta."; exit 1; }
    fi
fi

# Define required packages based on installation mode
if [ "$install_pyrosetta" = true ]; then
    required_packages=(pip pandas libgfortran5 matplotlib numpy biopython scipy pdbfixer openmm seaborn tqdm jupyter ffmpeg pyrosetta fsspec py3dmol chex dm-haiku dm-tree joblib ml-collections immutabledict optax jaxlib jax cuda-nvcc cudnn)
else
    required_packages=(pip pandas libgfortran5 matplotlib numpy biopython scipy pdbfixer openmm seaborn tqdm jupyter ffmpeg fsspec py3dmol chex dm-haiku dm-tree joblib ml-collections immutabledict optax jaxlib jax cuda-nvcc cudnn)
fi
missing_packages=()

# Check each package
for pkg in "${required_packages[@]}"; do
    conda list "$pkg" | grep -w "$pkg" >/dev/null 2>&1 || missing_packages+=("$pkg")
done

# If any packages are missing, output error and exit
if [ ${#missing_packages[@]} -ne 0 ]; then
    echo -e "Error: The following packages are missing from the environment:"
    for pkg in "${missing_packages[@]}"; do
        echo -e " - $pkg"
    done
    exit 1
fi

# install ColabDesign
echo -e "Installing ColabDesign\n"
pip3 install git+https://github.com/sokrypton/ColabDesign.git --no-deps || { echo -e "Error: Failed to install ColabDesign"; exit 1; }
python -c "import colabdesign" >/dev/null 2>&1 || { echo -e "Error: colabdesign module not found after installation"; exit 1; }

# install FreeSASA Python module
echo -e "Installing FreeSASA Python module\n"
pip3 install freesasa || { echo -e "Warning: Failed to install FreeSASA Python module via pip - FreeSASA SASA calculations will fall back to Biopython"; }
python -c "import freesasa" >/dev/null 2>&1 && echo -e "FreeSASA Python module installed successfully" || echo -e "Warning: FreeSASA Python module not available - using Biopython fallback for SASA"

# AlphaFold2 weights
echo -e "Downloading AlphaFold2 model weights \n"
params_dir="${install_dir}/params"
params_file="${params_dir}/alphafold_params_2022-12-06.tar"

# download AF2 weights
mkdir -p "${params_dir}" || { echo -e "Error: Failed to create weights directory"; exit 1; }
wget -O "${params_file}" "https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar" || { echo -e "Error: Failed to download AlphaFold2 weights"; exit 1; }
[ -s "${params_file}" ] || { echo -e "Error: Could not locate downloaded AlphaFold2 weights"; exit 1; }

# extract AF2 weights
tar tf "${params_file}" >/dev/null 2>&1 || { echo -e "Error: Corrupt AlphaFold2 weights download"; exit 1; }
tar -xvf "${params_file}" -C "${params_dir}" || { echo -e "Error: Failed to extract AlphaFold2weights"; exit 1; }
[ -f "${params_dir}/params_model_5_ptm.npz" ] || { echo -e "Error: Could not locate extracted AlphaFold2 weights"; exit 1; }
rm "${params_file}" || { echo -e "Warning: Failed to remove AlphaFold2 weights archive"; }

# chmod executables
echo -e "Changing permissions for executables\n"
chmod +x "${install_dir}/functions/dssp" || { echo -e "Error: Failed to chmod dssp"; exit 1; }

# chmod binaries in functions/ (sc, dssp, FASPR)
if [ -f "${install_dir}/functions/sc" ]; then
    chmod +x "${install_dir}/functions/sc" || { echo -e "Error: Failed to chmod sc"; exit 1; }
    echo -e "Made sc binary executable for shape complementarity calculations"
else
    echo -e "Warning: sc binary not found at ${install_dir}/functions/sc - shape complementarity will use placeholder values"
fi

if [ -f "${install_dir}/functions/dssp" ]; then
    chmod +x "${install_dir}/functions/dssp" || { echo -e "Error: Failed to chmod dssp"; exit 1; }
    echo -e "Made dssp binary executable"
else
    echo -e "Warning: dssp binary not found at ${install_dir}/functions/dssp - DSSP computations may be unavailable"
fi

if [ -f "${install_dir}/functions/FASPR" ]; then
    chmod +x "${install_dir}/functions/FASPR" || { echo -e "Error: Failed to chmod FASPR"; exit 1; }
    echo -e "Made FASPR binary executable"
else
    echo -e "Warning: FASPR binary not found at ${install_dir}/functions/FASPR - side-chain repacking will be unavailable"
fi

# Only setup DAlphaBall.gcc if installing PyRosetta
if [ "$install_pyrosetta" = true ]; then
    chmod +x "${install_dir}/functions/DAlphaBall.gcc" || { echo -e "Error: Failed to chmod DAlphaBall.gcc"; exit 1; }
else
    echo -e "Skipping DAlphaBall.gcc setup as PyRosetta is not being installed"
fi

# finish
conda deactivate
echo -e "BindCraft environment set up\n"

############################################################################################################
############################################################################################################
################## cleanup
echo -e "Cleaning up ${pkg_manager} temporary files to save space\n"
$pkg_manager clean -a -y
echo -e "$pkg_manager cleaned up\n"

################## finish script
t=$SECONDS 
echo -e "Successfully finished BindCraft installation!\n"
echo -e "Activate environment using command: \"$pkg_manager activate BindCraft\""
echo -e "\n"
echo -e "Installation took $(($t / 3600)) hours, $((($t / 60) % 60)) minutes and $(($t % 60)) seconds."