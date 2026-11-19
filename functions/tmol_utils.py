"""PyRosetta-free relaxation and interface scoring backed by tmol and fallbacks.

Uses ``tmol`` (executed in an isolated worker subprocess by default so PyTorch
does not collide with JAX/Metal/XLA contexts) for FastRelax relaxation,
``binder_score``, separated ``interface_dG``, ``interface_dG_SASA_ratio``, and
``interface_interface_hbonds`` / ``interface_hbond_percentage``. All remaining
metrics (``sc-rs`` shape complementarity, FreeSASA/Biopython SASA metrics,
interface residue accounting, and ``interface_packstat`` / ``interface_delta_unsat_hbonds``
placeholders) fall back to ``functions.pr_alternative_utils``.

Important: This module deliberately performs no top-level ``torch`` or ``tmol``
imports so importing ``functions`` in the JAX process remains lightweight and free
of PyTorch initialization.
"""

import atexit
import collections
import importlib.util
import json
import logging
import os
import pathlib
import selectors
import shutil
import string
import subprocess
import sys
import tempfile
import threading
import time

from Bio.PDB import PDBIO, PDBParser, Polypeptide

from . import pr_alternative_utils as alt
from .biopython_utils import (
    biopython_align_all_ca,
    compute_target_segment_lengths,
    merge_chains_into_single,
    split_chain_into_subchains,
)
from .generic_utils import clean_pdb
from .logging_utils import vprint

_TMOL_AVAILABILITY_CACHE = {}
_PERSISTENT_WORKERS = {}
_PERSISTENT_WORKERS_LOCK = threading.Lock()


class _PersistentTmolWorker:
    """Long-lived PyTorch/tmol worker subprocess communicating via JSON lines."""

    def __init__(self, python_exe):
        self.python_exe = os.path.abspath(python_exe)
        self._proc = None
        self._lock = threading.Lock()
        self._stderr_lines = collections.deque(maxlen=200)
        self._stderr_thread = None

    @property
    def pid(self):
        proc = self._proc
        if proc is not None and proc.poll() is None:
            return proc.pid
        return None

    def _drain_stderr(self, pipe):
        try:
            for raw_line in iter(pipe.readline, b""):
                text = raw_line.decode("utf-8", errors="replace")
                self._stderr_lines.append(text.rstrip())
                if logging.getLogger("functions").isEnabledFor(logging.DEBUG):
                    try:
                        sys.stderr.write(text)
                    except Exception:
                        pass
        except Exception:
            pass
        finally:
            try:
                pipe.close()
            except Exception:
                pass

    def _ensure_started(self):
        if self._proc is not None and self._proc.poll() is None:
            return
        self._cleanup_process()
        self._stderr_lines.clear()
        proc = subprocess.Popen(
            [self.python_exe, _tmol_worker_path(), "serve"],
            cwd=_repo_root(),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            bufsize=0,
        )
        self._proc = proc
        self._stderr_thread = threading.Thread(
            target=self._drain_stderr,
            args=(proc.stderr,),
            name=f"tmol-worker-stderr-{proc.pid}",
            daemon=True,
        )
        self._stderr_thread.start()

    def _cleanup_process(self, graceful=False):
        proc = self._proc
        self._proc = None
        if proc is None:
            return
        try:
            if graceful and proc.poll() is None and proc.stdin is not None:
                try:
                    proc.stdin.write(b'{"cmd": "shutdown"}\n')
                    proc.stdin.flush()
                except Exception:
                    pass
            if proc.stdin is not None:
                try:
                    proc.stdin.close()
                except Exception:
                    pass
            if proc.poll() is None:
                if graceful:
                    try:
                        proc.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        proc.terminate()
                        try:
                            proc.wait(timeout=3)
                        except subprocess.TimeoutExpired:
                            proc.kill()
                            proc.wait(timeout=3)
                else:
                    proc.terminate()
                    try:
                        proc.wait(timeout=2)
                    except subprocess.TimeoutExpired:
                        proc.kill()
                        proc.wait(timeout=2)
            else:
                proc.wait()
            for stream in (proc.stdout, proc.stderr):
                if stream is not None:
                    try:
                        stream.close()
                    except Exception:
                        pass
        except Exception:
            pass

    def _stderr_tail(self):
        return "\n".join(self._stderr_lines).strip()

    def _read_response_line(self, proc, timeout):
        stdout = proc.stdout
        if stdout is None:
            raise EOFError("tmol worker stdout closed")
        if timeout is None:
            line = stdout.readline()
            if not line:
                raise EOFError("tmol worker exited unexpectedly (EOF)")
            return line

        fd = stdout.fileno()
        buf = bytearray()
        deadline = time.monotonic() + float(timeout)
        with selectors.DefaultSelector() as sel:
            sel.register(fd, selectors.EVENT_READ)
            while True:
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    raise TimeoutError(f"tmol worker timed out after {timeout}s")
                events = sel.select(timeout=remaining)
                if not events:
                    raise TimeoutError(f"tmol worker timed out after {timeout}s")
                # Lockstep protocol sends one response per request; read 1 byte at a time
                # only for the final newline boundary if needed, or read small chunks?
                # Wait: reading 1-line JSON with os.read(fd, 1) avoids over-reading into
                # unbuffered FileIO and takes negligible time for ~250-byte JSON responses.
                ch = os.read(fd, 1)
                if not ch:
                    raise EOFError("tmol worker exited unexpectedly (EOF)")
                buf.extend(ch)
                if ch == b"\n":
                    return bytes(buf)

    def _send_once(self, payload, timeout):
        self._ensure_started()
        proc = self._proc
        req_bytes = (json.dumps(payload) + "\n").encode("utf-8")
        try:
            proc.stdin.write(req_bytes)
            proc.stdin.flush()
            raw_line = self._read_response_line(proc, timeout)
        except TimeoutError:
            self._cleanup_process(graceful=False)
            raise
        except (BrokenPipeError, EOFError, OSError) as exc:
            self._cleanup_process(graceful=False)
            raise exc

        resp = json.loads(raw_line.decode("utf-8"))
        if not isinstance(resp, dict):
            raise RuntimeError(f"Malformed tmol worker response: {resp!r}")
        return resp

    def request(self, payload, timeout=None):
        with self._lock:
            last_exc = None
            for attempt in range(2):
                try:
                    resp = self._send_once(payload, timeout)
                    break
                except (BrokenPipeError, EOFError, OSError) as exc:
                    last_exc = exc
                    if attempt == 0:
                        vprint(
                            f"[tmol-Worker] Restarting persistent worker after unexpected disconnect: {exc}"
                        )
                        continue
                    tail = self._stderr_tail()
                    detail = f"{exc}\n{tail}".strip() if tail else str(exc)
                    raise RuntimeError(f"tmol persistent worker failed: {detail}") from exc
            else:
                raise RuntimeError(f"tmol persistent worker failed: {last_exc}")

        if not resp.get("ok"):
            if logging.getLogger("functions").isEnabledFor(logging.DEBUG) and resp.get("traceback"):
                sys.stderr.write(resp["traceback"])
            raise RuntimeError(resp.get("error") or "tmol worker command failed")
        return resp["result"]

    def shutdown(self):
        with self._lock:
            self._cleanup_process(graceful=True)


def _get_persistent_worker(python_exe):
    key = os.path.abspath(python_exe)
    with _PERSISTENT_WORKERS_LOCK:
        worker = _PERSISTENT_WORKERS.get(key)
        if worker is None:
            worker = _PersistentTmolWorker(key)
            _PERSISTENT_WORKERS[key] = worker
        return worker


def get_tmol_worker_pid():
    """Return PID of the active persistent tmol worker for ``resolve_tmol_python()``, if any."""
    key = os.path.abspath(resolve_tmol_python())
    with _PERSISTENT_WORKERS_LOCK:
        worker = _PERSISTENT_WORKERS.get(key)
        return worker.pid if worker is not None else None


def shutdown_tmol_worker():
    """Gracefully stop any active persistent tmol worker subprocesses."""
    with _PERSISTENT_WORKERS_LOCK:
        workers = list(_PERSISTENT_WORKERS.values())
        _PERSISTENT_WORKERS.clear()
    for worker in workers:
        worker.shutdown()


atexit.register(shutdown_tmol_worker)


def _use_persistent_worker(persistent=None):
    if persistent is not None:
        return bool(persistent)
    env_val = os.environ.get("BINDCRAFT_TMOL_PERSISTENT", "1").strip().lower()
    return env_val not in ("0", "false", "off", "no")


def _payload_to_cli_args(command, payload):
    if command == "relax":
        args = [
            "relax",
            "--input",
            payload["input"],
            "--output",
            payload["output"],
            "--num-repeats",
            str(int(payload.get("num_repeats", 1))),
            "--constraint-weight",
            str(float(payload.get("constraint_weight", 1.0))),
        ]
        if not payload.get("constrain_coords", True):
            args.append("--no-constrain-coords")
        if payload.get("device") is not None:
            args.extend(["--device", str(payload["device"])])
        if payload.get("verbose"):
            args.append("--verbose")
        return args
    if command == "relax_batch":
        args = [
            "relax-batch",
            "--inputs",
            *[str(p) for p in payload["inputs"]],
            "--outputs",
            *[str(p) for p in payload["outputs"]],
            "--num-repeats",
            str(int(payload.get("num_repeats", 1))),
            "--constraint-weight",
            str(float(payload.get("constraint_weight", 1.0))),
        ]
        if not payload.get("constrain_coords", True):
            args.append("--no-constrain-coords")
        if payload.get("device") is not None:
            args.extend(["--device", str(payload["device"])])
        if payload.get("verbose"):
            args.append("--verbose")
        return args
    if command == "score":
        args = [
            "score",
            "--input",
            payload["input"],
            "--binder-chain",
            str(payload.get("binder_chain", "B")),
            "--target-chain",
            str(payload.get("target_chain", "A")),
            "--hbond-cutoff",
            str(float(payload.get("hbond_cutoff", -0.1))),
            "--seed",
            str(int(payload.get("seed", 0))),
        ]
        if not payload.get("pack_separated", True):
            args.append("--no-pack-separated")
        if payload.get("device") is not None:
            args.extend(["--device", str(payload["device"])])
        return args
    if command == "check":
        args = ["check"]
        if payload.get("pdb") is not None:
            args.extend(["--pdb", str(payload["pdb"])])
        if payload.get("device") is not None:
            args.extend(["--device", str(payload["device"])])
        return args
    raise ValueError(f"Unknown worker command: {command!r}")


def _invoke_tmol_worker(command, payload, timeout=None, persistent=None):
    """Run ``command`` in either a persistent daemon worker or a one-shot subprocess."""
    python_exe = resolve_tmol_python()
    if _use_persistent_worker(persistent):
        worker = _get_persistent_worker(python_exe)
        return worker.request({"cmd": command, **payload}, timeout=timeout)
    return _run_tmol_worker_subprocess(
        _payload_to_cli_args(command, payload),
        tmol_python=python_exe,
        timeout=timeout,
    )


def _is_executable_file(path):
    return bool(path) and os.path.isfile(path) and os.access(path, os.X_OK)


def resolve_tmol_python():
    """Return the Python interpreter path used to run tmol workers.

    Checks ``TMOL_PYTHON`` first, then ``$CONDA_PREFIX/tmol-venv/bin/python``
    (for Linux CUDA installs isolating PyTorch from JAX), and finally falls
    back to ``sys.executable``.
    """
    env_python = os.environ.get("TMOL_PYTHON", "").strip()
    if env_python and _is_executable_file(env_python):
        return env_python

    conda_prefix = os.environ.get("CONDA_PREFIX", "").strip()
    if conda_prefix:
        venv_python = os.path.join(conda_prefix, "tmol-venv", "bin", "python")
        if _is_executable_file(venv_python):
            return venv_python

    return sys.executable


def is_tmol_available(force_refresh=False):
    """Return True if ``tmol`` and ``torch`` are importable in ``resolve_tmol_python()``."""
    tmol_python = resolve_tmol_python()
    if not force_refresh and tmol_python in _TMOL_AVAILABILITY_CACHE:
        return _TMOL_AVAILABILITY_CACHE[tmol_python]

    available = False
    if os.path.abspath(tmol_python) == os.path.abspath(sys.executable):
        available = (
            importlib.util.find_spec("tmol") is not None
            and importlib.util.find_spec("torch") is not None
        )
    else:
        try:
            proc = subprocess.run(
                [
                    tmol_python,
                    "-c",
                    "import importlib.util, sys; "
                    "sys.exit(0 if (importlib.util.find_spec('tmol') "
                    "and importlib.util.find_spec('torch')) else 1)",
                ],
                capture_output=True,
                text=True,
                timeout=15,
            )
            available = proc.returncode == 0
        except Exception:
            available = False

    _TMOL_AVAILABILITY_CACHE[tmol_python] = available
    return available


def _tmol_worker_path():
    return str(pathlib.Path(__file__).resolve().with_name("tmol_worker.py"))


def _repo_root():
    return str(pathlib.Path(__file__).resolve().parents[1])


def _run_tmol_worker_subprocess(worker_args, tmol_python=None, timeout=None):
    """Execute ``functions/tmol_worker.py`` via ``resolve_tmol_python()`` and parse result JSON."""
    python_exe = tmol_python or resolve_tmol_python()
    worker = _tmol_worker_path()
    want_verbose = logging.getLogger("functions").isEnabledFor(logging.DEBUG)

    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as tf:
        result_json_path = tf.name

    try:
        cmd = [python_exe, worker, *worker_args, "--result-json", result_json_path]
        proc = subprocess.run(
            cmd,
            cwd=_repo_root(),
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        if want_verbose:
            if proc.stdout:
                sys.stdout.write(proc.stdout)
            if proc.stderr:
                sys.stderr.write(proc.stderr)
        if proc.returncode != 0:
            err_detail = (proc.stderr or proc.stdout or "").strip()
            raise RuntimeError(
                f"tmol_worker exited with code {proc.returncode}: {err_detail}"
            )
        with open(result_json_path, "r", encoding="utf-8") as handle:
            return json.load(handle)
    finally:
        try:
            if os.path.isfile(result_json_path):
                os.remove(result_json_path)
        except Exception:
            pass


def _capture_residue_bfactors(pdb_file_path):
    """Record original per-residue B-factors (CA or first atom) keyed by (chain, residue_id)."""
    bfactors = {}
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("original", pdb_file_path)
        for model in structure:
            for chain in model:
                for residue in chain:
                    if not Polypeptide.is_aa(residue, standard=True):
                        continue
                    atom = residue["CA"] if "CA" in residue else next(residue.get_atoms(), None)
                    if atom is not None:
                        bfactors[(chain.id, residue.id)] = float(atom.get_bfactor())
    except Exception:
        return {}
    return bfactors


def _restore_residue_bfactors(pdb_file_path, bfactors):
    """Apply saved residue B-factors onto all atoms of matching residues in ``pdb_file_path``."""
    if not bfactors:
        return
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("relaxed", pdb_file_path)
        changed = False
        for model in structure:
            for chain in model:
                for residue in chain:
                    val = bfactors.get((chain.id, residue.id))
                    if val is not None:
                        for atom in residue:
                            atom.set_bfactor(val)
                        changed = True
        if changed:
            io = PDBIO()
            io.set_structure(structure)
            io.save(pdb_file_path)
    except Exception:
        pass


def _prepare_deconcat_pdb(pdb_file_path):
    """Split concatenated chain 'A' into target segments if BINDCRAFT env context requires it."""
    starting_pdb = os.environ.get("BINDCRAFT_STARTING_PDB", "").strip()
    target_chains = os.environ.get("BINDCRAFT_TARGET_CHAINS", "").strip()
    if not (starting_pdb and os.path.isfile(starting_pdb) and target_chains):
        return pdb_file_path, None, None

    try:
        lengths = compute_target_segment_lengths(starting_pdb, target_chains)
        if not lengths or sum(1 for length in lengths if length > 0) < 2:
            return pdb_file_path, None, None

        pool = [c for c in string.ascii_uppercase if c not in ("A", "B")]
        pool += list("0123456789")
        pool += list(string.ascii_lowercase)
        n_segments = len(lengths)
        if len(pool) < n_segments:
            raise RuntimeError(f"Insufficient temporary chain IDs for {n_segments} target segments")
        new_chain_ids = pool[:n_segments]

        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tf:
            deconcat_tmp = tf.name

        split_chain_into_subchains(
            pdb_file_path,
            source_chain_id="A",
            subchain_lengths=lengths,
            new_chain_ids=new_chain_ids,
            output_path=deconcat_tmp,
        )
        if not os.path.isfile(deconcat_tmp) or os.path.getsize(deconcat_tmp) == 0:
            try:
                if os.path.isfile(deconcat_tmp):
                    os.remove(deconcat_tmp)
            except Exception:
                pass
            return pdb_file_path, None, None

        return deconcat_tmp, deconcat_tmp, (new_chain_ids, "A")
    except Exception as exc:
        vprint(f"[tmol-Relax] De-concatenation skipped: {exc}")
        return pdb_file_path, None, None


def _prepare_single_relax_item(src_pdb, dst_pdb):
    """Prepare original B-factors, optional de-concatenation, and temp output path."""
    os.makedirs(os.path.dirname(dst_pdb) or ".", exist_ok=True)
    original_bfactors = _capture_residue_bfactors(src_pdb)
    work_in, deconcat_tmp, reconcat_spec = _prepare_deconcat_pdb(src_pdb)
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tf_out:
        work_out = tf_out.name
    return {
        "src_pdb": src_pdb,
        "dst_pdb": dst_pdb,
        "basename": os.path.basename(src_pdb),
        "original_bfactors": original_bfactors,
        "work_in": work_in,
        "deconcat_tmp": deconcat_tmp,
        "reconcat_spec": reconcat_spec,
        "work_out": work_out,
    }


def _finalize_single_relax_item(item):
    """Validate, re-concatenate, align, restore B-factors, clean, and move relaxed PDB."""
    work_out = item["work_out"]
    if not work_out or not os.path.isfile(work_out) or os.path.getsize(work_out) == 0:
        raise RuntimeError(f"tmol relax produced empty output for {item['basename']}")

    reconcat_spec = item["reconcat_spec"]
    if reconcat_spec is not None:
        src_chains, dest_chain = reconcat_spec
        vprint(f"[tmol-Relax] Re-concatenating segments {src_chains} back into chain {dest_chain}")
        merge_chains_into_single(
            work_out, src_chains, dest_chain_id=dest_chain, output_path=work_out
        )

    biopython_align_all_ca(item["src_pdb"], work_out)
    _restore_residue_bfactors(work_out, item["original_bfactors"])
    clean_pdb(work_out)
    shutil.move(work_out, item["dst_pdb"])
    item["work_out"] = None


def _cleanup_relax_item(item):
    for tmp_file in (item.get("deconcat_tmp"), item.get("work_out")):
        try:
            if tmp_file and os.path.isfile(tmp_file):
                os.remove(tmp_file)
        except Exception:
            pass


def tmol_relax(
    pdb_file_path,
    output_pdb_path,
    use_subprocess=True,
    fallback_to_openmm=True,
    num_repeats=1,
    constrain_coords=True,
    constraint_weight=1.0,
    timeout=None,
    device=None,
    perf_report=None,
    persistent=None,
):
    """Relax a PDB structure with tmol FastRelax, aligning and preserving B-factors.

    Falls back to ``openmm_relax_subprocess`` / ``openmm_relax`` on failure when
    ``fallback_to_openmm=True``.
    """
    t0 = time.time()
    basename = os.path.basename(pdb_file_path)
    src_pdb = os.path.abspath(pdb_file_path)
    dst_pdb = os.path.abspath(output_pdb_path)
    vprint(f"[tmol-Relax] Initiating tmol FastRelax for {basename} (subprocess={use_subprocess})")

    item = None
    try:
        item = _prepare_single_relax_item(src_pdb, dst_pdb)

        if use_subprocess:
            payload = {
                "input": item["work_in"],
                "output": item["work_out"],
                "num_repeats": int(num_repeats),
                "constrain_coords": bool(constrain_coords),
                "constraint_weight": float(constraint_weight),
                "device": device,
                "verbose": logging.getLogger("functions").isEnabledFor(logging.DEBUG),
            }
            meta = _invoke_tmol_worker("relax", payload, timeout=timeout, persistent=persistent)
        else:
            from . import tmol_worker

            meta = tmol_worker.run_tmol_relax(
                item["work_in"],
                item["work_out"],
                num_repeats=num_repeats,
                constrain_coords=constrain_coords,
                constraint_weight=constraint_weight,
                device=device,
                verbose=logging.getLogger("functions").isEnabledFor(logging.DEBUG),
            )

        _finalize_single_relax_item(item)

        elapsed = time.time() - t0
        if isinstance(perf_report, dict):
            perf_report.update(meta)
            perf_report["total_seconds"] = elapsed
            perf_report["used_fallback"] = False
        vprint(
            f"[tmol-Relax] Completed {basename} in {elapsed:.2f}s "
            f"(device={meta.get('device')}, mps_offload={meta.get('mps_offload')}, "
            f"E={meta.get('initial_energy', 0):.1f}->{meta.get('final_energy', 0):.1f})"
        )
        return meta
    except Exception as exc:
        if not fallback_to_openmm:
            raise
        print(f"[tmol-Relax] WARN: tmol FastRelax failed for {basename} ({exc}); falling back to OpenMM relax.")
        if isinstance(perf_report, dict):
            perf_report["used_fallback"] = True
            perf_report["fallback_reason"] = str(exc)
        if use_subprocess:
            return alt.openmm_relax_subprocess(src_pdb, dst_pdb, use_gpu_relax=True)
        return alt.openmm_relax(src_pdb, dst_pdb, use_gpu_relax=True)
    finally:
        if item is not None:
            _cleanup_relax_item(item)


def tmol_relax_subprocess(pdb_file_path, output_pdb_path, **kwargs):
    """Run ``tmol_relax`` inside a fresh isolated worker subprocess."""
    kwargs.pop("use_subprocess", None)
    return tmol_relax(pdb_file_path, output_pdb_path, use_subprocess=True, **kwargs)


def tmol_relax_batch(
    pdb_pairs,
    use_subprocess=True,
    fallback_to_openmm=True,
    num_repeats=1,
    constrain_coords=True,
    constraint_weight=1.0,
    timeout=None,
    device=None,
    perf_report=None,
    persistent=None,
):
    """Relax multiple ``(input_pdb, output_pdb)`` structures together in one PoseStack.

    Falls back cleanly to sequential ``tmol_relax`` (which preserves per-structure
    OpenMM fallback if needed) when batch execution fails.
    """
    pairs = [(os.path.abspath(src), os.path.abspath(dst)) for src, dst in pdb_pairs]
    if not pairs:
        return {}
    if len(pairs) == 1:
        return tmol_relax(
            pairs[0][0],
            pairs[0][1],
            use_subprocess=use_subprocess,
            fallback_to_openmm=fallback_to_openmm,
            num_repeats=num_repeats,
            constrain_coords=constrain_coords,
            constraint_weight=constraint_weight,
            timeout=timeout,
            device=device,
            perf_report=perf_report,
            persistent=persistent,
        )

    t0 = time.time()
    vprint(
        f"[tmol-RelaxBatch] Initiating batched tmol FastRelax for {len(pairs)} structures "
        f"(subprocess={use_subprocess})"
    )
    items = []
    try:
        for src, dst in pairs:
            items.append(_prepare_single_relax_item(src, dst))
        inputs = [item["work_in"] for item in items]
        outputs = [item["work_out"] for item in items]

        if use_subprocess:
            payload = {
                "inputs": inputs,
                "outputs": outputs,
                "num_repeats": int(num_repeats),
                "constrain_coords": bool(constrain_coords),
                "constraint_weight": float(constraint_weight),
                "device": device,
                "verbose": logging.getLogger("functions").isEnabledFor(logging.DEBUG),
            }
            meta = _invoke_tmol_worker("relax_batch", payload, timeout=timeout, persistent=persistent)
        else:
            from . import tmol_worker

            meta = tmol_worker.run_tmol_relax_batch(
                inputs,
                outputs,
                num_repeats=num_repeats,
                constrain_coords=constrain_coords,
                constraint_weight=constraint_weight,
                device=device,
                verbose=logging.getLogger("functions").isEnabledFor(logging.DEBUG),
            )

        for item in items:
            _finalize_single_relax_item(item)

        elapsed = time.time() - t0
        if isinstance(perf_report, dict):
            perf_report.update(meta)
            perf_report["total_seconds"] = elapsed
            perf_report["used_fallback"] = False
            perf_report["used_batch_fallback"] = False
        vprint(
            f"[tmol-RelaxBatch] Completed {len(pairs)} structures in {elapsed:.2f}s "
            f"(device={meta.get('device')}, mps_offload={meta.get('mps_offload')})"
        )
        return meta
    except Exception as exc:
        for item in items:
            _cleanup_relax_item(item)
        items = []
        if not fallback_to_openmm:
            raise
        print(
            f"[tmol-RelaxBatch] WARN: batched tmol FastRelax failed ({exc}); "
            "falling back to per-structure relaxation."
        )
        if isinstance(perf_report, dict):
            perf_report["used_batch_fallback"] = True
            perf_report["batch_fallback_reason"] = str(exc)
        results = []
        for src, dst in pairs:
            results.append(
                tmol_relax(
                    src,
                    dst,
                    use_subprocess=use_subprocess,
                    fallback_to_openmm=fallback_to_openmm,
                    num_repeats=num_repeats,
                    constrain_coords=constrain_coords,
                    constraint_weight=constraint_weight,
                    timeout=timeout,
                    device=device,
                    persistent=persistent,
                )
            )
        return {"fallback_results": results}
    finally:
        for item in items:
            _cleanup_relax_item(item)


def tmol_relax_batch_subprocess(pdb_pairs, **kwargs):
    """Run ``tmol_relax_batch`` inside an isolated worker subprocess."""
    kwargs.pop("use_subprocess", None)
    return tmol_relax_batch(pdb_pairs, use_subprocess=True, **kwargs)


def tmol_score_interface(
    pdb_file,
    binder_chain="B",
    target_chain="A",
    sasa_engine="auto",
    use_subprocess=True,
    fallback_to_alt=True,
    pack_separated=True,
    hbond_energy_cutoff=-0.1,
    timeout=None,
    device=None,
    seed=0,
    persistent=None,
):
    """Score interface using tmol physics where supported and pr_alternative_utils elsewhere.

    Preserves ``interface_packstat`` (0.65), ``interface_delta_unsat_hbonds`` (1),
    and ``interface_delta_unsat_hbonds_percentage`` (0.0) placeholders from
    ``pr_alternative_score_interface``, along with ``sc-rs`` shape complementarity,
    FreeSASA/Biopython SASA metrics, and Biopython interface residue lists.
    """
    t0 = time.time()
    basename = os.path.basename(pdb_file)
    vprint(f"[tmol-Score] Scoring {basename} (binder={binder_chain}, subprocess={use_subprocess})")

    interface_scores, interface_AA, interface_residues_pdb_ids_str = (
        alt.pr_alternative_score_interface(
            pdb_file,
            binder_chain=binder_chain,
            target_chain=target_chain,
            sasa_engine=sasa_engine,
        )
    )

    try:
        if use_subprocess:
            payload = {
                "input": os.path.abspath(pdb_file),
                "binder_chain": str(binder_chain),
                "target_chain": str(target_chain),
                "pack_separated": bool(pack_separated),
                "hbond_cutoff": float(hbond_energy_cutoff),
                "device": device,
                "seed": int(seed),
            }
            tmol_res = _invoke_tmol_worker("score", payload, timeout=timeout, persistent=persistent)
        else:
            from . import tmol_worker

            tmol_res = tmol_worker.run_tmol_score(
                pdb_file,
                binder_chain=binder_chain,
                target_chain=target_chain,
                pack_separated=pack_separated,
                hbond_energy_cutoff=hbond_energy_cutoff,
                device=device,
                seed=seed,
            )

        interface_dSASA = float(interface_scores.get("interface_dSASA") or 0.0)
        interface_nres = int(interface_scores.get("interface_nres") or 0)

        binder_score = float(tmol_res["binder_score"])
        interface_dG = float(tmol_res["interface_dG"])
        interface_hbonds = int(tmol_res["interface_interface_hbonds"])

        dG_sasa_ratio = (interface_dG / interface_dSASA) * 100.0 if interface_dSASA > 0.0 else 0.0
        hbond_pct = (interface_hbonds / interface_nres) * 100.0 if interface_nres > 0 else 0.0

        interface_scores["binder_score"] = round(binder_score, 2) + 0.0
        interface_scores["interface_dG"] = round(interface_dG, 2) + 0.0
        interface_scores["interface_dG_SASA_ratio"] = round(dG_sasa_ratio, 2) + 0.0
        interface_scores["interface_interface_hbonds"] = interface_hbonds
        interface_scores["interface_hbond_percentage"] = round(hbond_pct, 2) + 0.0

        vprint(
            f"[tmol-Score] Completed {basename} in {time.time()-t0:.2f}s "
            f"(binder_score={interface_scores['binder_score']}, dG={interface_scores['interface_dG']}, "
            f"hbonds={interface_hbonds}, mps_offload={tmol_res.get('mps_offload')})"
        )
        return interface_scores, interface_AA, interface_residues_pdb_ids_str
    except Exception as exc:
        if not fallback_to_alt:
            raise
        print(f"[tmol-Score] WARN: tmol scoring failed for {basename} ({exc}); falling back to alternative placeholders.")
        return interface_scores, interface_AA, interface_residues_pdb_ids_str


def tmol_score_interface_subprocess(pdb_file, binder_chain="B", **kwargs):
    """Run ``tmol_score_interface`` using a fresh isolated worker subprocess."""
    kwargs.pop("use_subprocess", None)
    return tmol_score_interface(pdb_file, binder_chain=binder_chain, use_subprocess=True, **kwargs)
