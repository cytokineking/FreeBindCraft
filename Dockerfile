FROM nvidia/cuda:12.1.1-cudnn8-runtime-ubuntu22.04

LABEL org.opencontainers.image.source="https://github.com/cytokineking/FreeBindCraft"
LABEL org.opencontainers.image.description="FreeBindCraft GPU (no PyRosetta)"
LABEL org.opencontainers.image.licenses="MIT"

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Etc/UTC

# OS dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      bash \
      ca-certificates \
      curl \
      git \
      rsync \
      libgfortran5 \
      tmux \
      wget \
      build-essential \
      pkg-config \
      procps \
      unzip && \
    rm -rf /var/lib/apt/lists/*

# Install OpenCL ICD loader and tools; register NVIDIA OpenCL ICD
RUN apt-get update && apt-get install -y --no-install-recommends \
    ocl-icd-libopencl1 clinfo && \
    rm -rf /var/lib/apt/lists/*
RUN mkdir -p /etc/OpenCL/vendors && \
    echo "libnvidia-opencl.so.1" > /etc/OpenCL/vendors/nvidia.icd

# Install Miniforge (Conda) at /miniforge3
ENV CONDA_DIR=/miniforge3
RUN wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p ${CONDA_DIR} && \
    rm -f /tmp/miniforge.sh

# Put conda on PATH
ENV PATH=${CONDA_DIR}/bin:${PATH}

# Improve conda robustness and cleanup
RUN conda config --set channel_priority strict && \
    conda config --set always_yes yes && \
    conda update -n base -c conda-forge conda && \
    conda clean -afy

# Create workdir and copy project
WORKDIR /app
COPY . /app

# Ensure helper binaries are executable (also handled by installer)
RUN chmod +x /app/functions/dssp || true && \
    chmod +x /app/functions/sc || true

# Build environment and download AF2 weights without PyRosetta
# Match CUDA to base image; installer pins jax/jaxlib=0.6.0
RUN bash -lc 'source ${CONDA_DIR}/etc/profile.d/conda.sh && \
    bash /app/install_bindcraft.sh --pkg_manager conda --cuda 12.1 --no-pyrosetta'

# Default environment
ENV PATH=${CONDA_DIR}/envs/BindCraft/bin:${CONDA_DIR}/bin:${PATH} \
    LD_LIBRARY_PATH=${CONDA_DIR}/envs/BindCraft/lib:${LD_LIBRARY_PATH} \
    PYTHONUNBUFFERED=1 \
    BINDCRAFT_HOME=/app

# Prefer OpenCL (fallback to CUDA) in OpenMM by default
ENV OPENMM_PLATFORM_ORDER=OpenCL,CUDA \
    OPENMM_DEFAULT_PLATFORM=OpenCL

# Modal-compatible entrypoint that execs args
COPY docker-entrypoint.sh /usr/local/bin/bindcraft-entrypoint.sh
RUN chmod +x /usr/local/bin/bindcraft-entrypoint.sh
ENTRYPOINT ["/usr/local/bin/bindcraft-entrypoint.sh"]

# Default command prints help
CMD ["python", "bindcraft.py", "--help"]
