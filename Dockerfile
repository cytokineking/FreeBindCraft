# The pinned pip profile supplies CUDA libraries and cuDNN.
FROM nvidia/cuda:12.9.1-base-ubuntu22.04

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
# The installer supplies the pinned pip CUDA 12 profile; builds need no GPU.
# Allow toggling PyRosetta install at build-time
ARG WITH_PYROSETTA=false
ENV WITH_PYROSETTA=${WITH_PYROSETTA}
RUN bash -lc 'source ${CONDA_DIR}/etc/profile.d/conda.sh && \
    EXTRA=""; if [ "${WITH_PYROSETTA}" != "true" ]; then EXTRA="--no-pyrosetta"; fi; \
    bash /app/install_bindcraft.sh --pkg_manager conda --skip-gpu-check ${EXTRA} && \
    conda clean -afy && rm -rf /root/.cache/pip'

# Default environment
ENV PATH=${CONDA_DIR}/envs/BindCraft/bin:${CONDA_DIR}/bin:${PATH} \
    LD_LIBRARY_PATH="" \
    JAX_PLATFORMS=cuda \
    PYTHONUNBUFFERED=1 \
    BINDCRAFT_HOME=/app

# Use the CUDA plugin shipped with the tested OpenMM profile.
ENV OPENMM_PLATFORM_ORDER=CUDA \
    OPENMM_DEFAULT_PLATFORM=CUDA

# Modal-compatible entrypoint that execs args
COPY docker-entrypoint.sh /usr/local/bin/bindcraft-entrypoint.sh
RUN chmod +x /usr/local/bin/bindcraft-entrypoint.sh
ENTRYPOINT ["/usr/local/bin/bindcraft-entrypoint.sh"]

# Default command prints help
CMD ["python", "bindcraft.py", "--help"]
