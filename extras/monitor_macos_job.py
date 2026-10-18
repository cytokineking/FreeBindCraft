#!/usr/bin/env python3
"""Run a real job and record macOS physical footprint, RSS and system memory.

Physical footprint and RSS do not separately measure Metal allocations. System
pressure, compression, and swap include other applications on the Mac.
"""
import argparse
from collections import Counter
import ctypes
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import re
import signal
import subprocess
import time

import psutil


class RusageInfoV4(ctypes.Structure):
    # Public macOS SDK sys/resource.h, RUSAGE_INFO_V4.
    _fields_ = [("ri_uuid", ctypes.c_uint8 * 16)] + [
        ("ri_" + name, ctypes.c_uint64) for name in (
            "user_time system_time pkg_idle_wkups interrupt_wkups pageins "
            "wired_size resident_size phys_footprint proc_start_abstime proc_exit_abstime "
            "child_user_time child_system_time child_pkg_idle_wkups child_interrupt_wkups "
            "child_pageins child_elapsed_abstime diskio_bytesread diskio_byteswritten "
            "cpu_time_qos_default cpu_time_qos_maintenance cpu_time_qos_background "
            "cpu_time_qos_utility cpu_time_qos_legacy cpu_time_qos_user_initiated "
            "cpu_time_qos_user_interactive billed_system_time serviced_system_time "
            "logical_writes lifetime_max_phys_footprint instructions cycles billed_energy "
            "serviced_energy interval_max_phys_footprint runnable_time"
        ).split()
    ]


def read_footprint(pid):
    libproc = ctypes.CDLL("/usr/lib/libproc.dylib", use_errno=True)
    libproc.proc_pid_rusage.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
    libproc.proc_pid_rusage.restype = ctypes.c_int
    usage = RusageInfoV4()
    if libproc.proc_pid_rusage(pid, 4, ctypes.byref(usage)) != 0:
        return {}
    return {"physical_footprint_bytes": usage.ri_phys_footprint,
            "lifetime_peak_physical_footprint_bytes": usage.ri_lifetime_max_phys_footprint}


def read_system():
    vm = psutil.virtual_memory()
    swap = psutil.swap_memory()
    result = {"available_bytes": vm.available, "used_bytes": vm.used,
              "wired_bytes": getattr(vm, "wired", None), "swap_used_bytes": swap.used}
    pressure = subprocess.run(["sysctl", "-n", "kern.memorystatus_vm_pressure_level"],
                              capture_output=True, text=True)
    result["pressure_level"] = int(pressure.stdout.strip()) if pressure.returncode == 0 else None
    stats = subprocess.run(["vm_stat"], capture_output=True, text=True, check=True).stdout
    page_size = int(re.search(r"page size of (\d+) bytes", stats).group(1))
    for label, key in (("Pages occupied by compressor", "compressed_bytes"),
                       ("Pages stored in compressor", "uncompressed_bytes_in_compressor"),
                       ("Swapouts", "swapouts_bytes")):
        match = re.search(re.escape(label) + r":\s+(\d+)", stats)
        result[key] = int(match.group(1)) * page_size if match else None
    return result


def read_tree(pid):
    rows = []
    try:
        parent = psutil.Process(pid)
        processes = [parent, *parent.children(recursive=True)]
    except psutil.NoSuchProcess:
        processes = []
    for process in processes:
        try:
            rows.append({"pid": process.pid, "name": process.name(),
                         "rss_bytes": process.memory_info().rss,
                         **read_footprint(process.pid)})
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--interval", type=float, default=1.0)
    parser.add_argument("--max-tree-rss-gib", type=float)
    parser.add_argument("--max-hours", type=float, default=3)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command or args.interval <= 0:
        parser.error("A command and positive sampling interval are required")
    if platform.system() != "Darwin":
        parser.error("This monitor requires macOS")
    destination = args.output_dir.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    if (destination / "memory.jsonl").exists():
        parser.error("Output directory already contains a run; choose a new directory")
    baseline = read_system()
    metadata = {"command": command, "cwd": str(Path.cwd()),
                "started_utc": datetime.now(timezone.utc).isoformat(),
                "macos": platform.mac_ver()[0], "machine": platform.machine(),
                "physical_memory_bytes": psutil.virtual_memory().total,
                "interval_seconds": args.interval, "baseline": baseline,
                "rss_note": "Sampled sum of process-tree RSS; shared pages may be counted twice. Metal allocations are not isolated.",
                "footprint_note": "macOS proc_pid_rusage physical footprint and per-process lifetime high-water mark; not a separate GPU allocation measurement.",
                "system_note": "System memory, compression, and swap include other running applications."}
    (destination / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    started = time.monotonic()
    samples = []
    stop_reason = None
    critical_since = None
    with (destination / "job.log").open("w") as log, (destination / "memory.jsonl").open("w") as stream:
        child = subprocess.Popen(["/usr/bin/time", "-l", *command], stdout=log,
                                 stderr=subprocess.STDOUT, start_new_session=True)
        (destination / "pid").write_text(str(child.pid) + "\n")
        while True:
            elapsed = time.monotonic() - started
            system = read_system()
            tree = read_tree(child.pid)
            sample = {"elapsed_seconds": round(elapsed, 3), "system": system,
                      "tree_rss_bytes": sum(row["rss_bytes"] for row in tree),
                      "tree_physical_footprint_bytes": sum(row.get("physical_footprint_bytes", 0) for row in tree),
                      "processes": tree}
            samples.append(sample)
            stream.write(json.dumps(sample) + "\n")
            stream.flush()
            if len(samples) % 10 == 1:
                (destination / "latest.json").write_text(json.dumps(sample, indent=2) + "\n")
            if child.poll() is not None:
                break
            if system["pressure_level"] == 4:
                critical_since = critical_since or time.monotonic()
            else:
                critical_since = None
            if args.max_tree_rss_gib and sample["tree_rss_bytes"] > args.max_tree_rss_gib * 1024**3:
                stop_reason = "process-tree RSS limit"
            elif critical_since and time.monotonic() - critical_since > 60 and system["swap_used_bytes"] - baseline["swap_used_bytes"] > 4 * 1024**3:
                stop_reason = "critical memory pressure for 60 seconds with over 4 GiB additional swap"
            elif elapsed > args.max_hours * 3600:
                stop_reason = "elapsed time limit"
            if stop_reason:
                os.killpg(child.pid, signal.SIGTERM)
                try:
                    child.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid, signal.SIGKILL)
                    child.wait()
                break
            time.sleep(args.interval)
    final = read_system()
    summary = {**metadata, "returncode": child.returncode, "stop_reason": stop_reason,
               "elapsed_seconds": time.monotonic() - started, "samples": len(samples),
               "peak_tree_rss_bytes": max(row["tree_rss_bytes"] for row in samples),
               "peak_sampled_tree_physical_footprint_bytes": max(row["tree_physical_footprint_bytes"] for row in samples),
               "maximum_process_lifetime_physical_footprint_bytes": max(
                   (process.get("lifetime_peak_physical_footprint_bytes", 0)
                    for row in samples for process in row["processes"]), default=0),
               "minimum_system_available_bytes": min(row["system"]["available_bytes"] for row in samples),
               "peak_system_compressed_bytes": max(row["system"]["compressed_bytes"] for row in samples),
               "peak_system_swap_used_bytes": max(row["system"]["swap_used_bytes"] for row in samples),
               "pressure_samples": dict(Counter(str(row["system"]["pressure_level"]) for row in samples)),
               "final_system": final}
    summary["peak_swap_growth_bytes"] = max(0, summary["peak_system_swap_used_bytes"] - baseline["swap_used_bytes"])
    summary["final_swap_change_bytes"] = final["swap_used_bytes"] - baseline["swap_used_bytes"]
    log_text = (destination / "job.log").read_text(errors="replace")
    # /usr/bin/time's "peak memory footprint" can describe the time wrapper
    # itself; use proc_pid_rusage on the actual job process instead.
    for label, key in (("maximum resident set size", "time_maximum_rss_bytes"),):
        match = re.search(r"^\s*(\d+)\s+" + re.escape(label) + r"\s*$", log_text, re.MULTILINE)
        summary[key] = int(match.group(1)) if match else None
    (destination / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (destination / "exit_status").write_text(str(child.returncode) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    raise SystemExit(child.returncode if child.returncode >= 0 else 128 - child.returncode)


if __name__ == "__main__":
    main()
