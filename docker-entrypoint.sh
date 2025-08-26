#!/usr/bin/env bash
set -euo pipefail

# Increase file descriptor limits if possible (BindCraft can open many files)
if command -v prlimit >/dev/null 2>&1; then
  prlimit --pid $$ --nofile=65536:65536 || true
else
  # Fallback: try ulimit within shell
  ulimit -n 65536 || true
fi

# Ensure /app is working directory
cd /app

# Exec passed command (required for Modal ENTRYPOINT compatibility)
exec "$@"
