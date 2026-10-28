"""Metal execution selection and an observable, conservative memory preflight.

Admission combines a necessary score-buffer bound and an empirical working-set
estimate. The native runtime remains responsible for actual allocations and
pressure handling; the estimate is not a guarantee.
"""
import os

_GIB = 1024 ** 3


def select_global_policy(global_config, total_residues, backend):
    """Apply the global policy before creating the model's final callables."""
    total_residues = int(total_residues)
    if total_residues <= 0:
        raise ValueError('The model must contain at least one residue')
    chunking = backend != 'metal' or os.environ.get('FREEBINDCRAFT_METAL_CHUNKING') == '1'
    desired = 4 if chunking and total_residues > 384 else None
    changed = global_config.subbatch_size != desired
    global_config.subbatch_size = desired
    return changed, {
        'backend': backend,
        'total_residues': total_residues,
        'global_subbatch_size': desired,
    }


def assess_metal_memory(total_residues, num_heads, element_bytes, *, total_bytes,
                        available_bytes, process_bytes=None, device_limit_bytes=None,
                        device_used_bytes=None, system_claimed_bytes=None,
                        system_ceiling_bytes=None, process_ceiling_bytes=None,
                        retained_model_bytes=None, nonreusable_process_bytes=None,
                        reclaimable_process_cache_bytes=None,
                        reclaimable_claimed_cache_bytes=None,
                        reported_device_cache_bytes=None,
                        rematerialization=True, recycles=None,
                        recycle_mode='last', model_kind=None):
    """Estimate incremental headroom without counting retained memory twice.

    Process bytes are macOS physical footprint, not RSS. System claimed bytes
    use the native governor's wired + internal - purgeable + compressed metric.
    Retained and reclaimable arguments are *physical bytes already included* in
    the corresponding measurement, not logical tensor sizes or cache estimates.
    Unclassified process memory is never credited as reusable: if its unknown
    ownership is the only reason a limit might be exceeded, the answer is
    unknown and the live runtime governor makes the actual allocation decision.
    """
    total_residues, num_heads, element_bytes = map(int, (total_residues, num_heads, element_bytes))
    if min(total_residues, num_heads, element_bytes) <= 0:
        raise ValueError('Positive residue, head and element sizes are required')

    def optional_bytes(value):
        return max(0, int(value)) if value is not None else None

    total_bytes, available_bytes, process_bytes = map(optional_bytes, (total_bytes, available_bytes, process_bytes))
    device_limit_bytes, device_used_bytes = map(optional_bytes, (device_limit_bytes, device_used_bytes))
    system_claimed_bytes, system_ceiling_bytes, process_ceiling_bytes = map(
        optional_bytes, (system_claimed_bytes, system_ceiling_bytes, process_ceiling_bytes))
    retained_model_bytes = optional_bytes(retained_model_bytes)
    nonreusable_process_bytes = optional_bytes(nonreusable_process_bytes)
    reclaimable_process_cache_bytes = optional_bytes(reclaimable_process_cache_bytes)
    reclaimable_claimed_cache_bytes = optional_bytes(reclaimable_claimed_cache_bytes)
    reported_device_cache_bytes = optional_bytes(reported_device_cache_bytes)
    recycles = int(recycles) if recycles is not None else None
    score_bytes = total_residues ** 3 * num_heads * element_bytes
    ratio = total_residues / 406.0
    scale = max(ratio ** 2, ratio ** 3)
    dtype_heads = max(1.0, (element_bytes / 2.0) * (num_heads / 4.0))
    calibrated = bool(rematerialization) and recycle_mode == 'last'
    estimate = int((6 + 24 * scale * dtype_heads) * _GIB) if calibrated else None

    # Do not invent ownership or exceed the footprint measurement. Only the
    # supplied, physically verified subset is credited as retained model data.
    retained = min(process_bytes, retained_model_bytes or 0) if process_bytes is not None else 0
    remaining = max(0, process_bytes - retained) if process_bytes is not None else None
    nonreusable = min(remaining, nonreusable_process_bytes or 0) if remaining is not None else 0
    remaining_after_live = max(0, remaining - nonreusable) if remaining is not None else None
    cache = min(remaining_after_live, reclaimable_process_cache_bytes or 0) if remaining_after_live is not None else 0
    unclassified = remaining_after_live - cache if remaining_after_live is not None else None
    claimed_cache = min(cache, reclaimable_claimed_cache_bytes or 0,
                        system_claimed_bytes or 0)
    reserve = max(2 * _GIB, min(8 * _GIB, total_bytes // 10)) if total_bytes else None
    result = {
        'status': 'unknown', 'score_tensor_bytes': score_bytes,
        'total_bytes': total_bytes, 'available_bytes': available_bytes,
        'process_footprint_bytes': process_bytes, 'system_claimed_bytes': system_claimed_bytes,
        'retained_model_footprint_bytes': retained_model_bytes,
        'credited_retained_model_bytes': retained,
        'known_nonreusable_process_bytes': nonreusable_process_bytes,
        'reclaimable_process_cache_bytes': reclaimable_process_cache_bytes,
        'reclaimable_claimed_cache_bytes': reclaimable_claimed_cache_bytes,
        'unclassified_process_footprint_bytes': unclassified,
        'reported_device_cache_bytes': reported_device_cache_bytes,
        'device_limit_bytes': device_limit_bytes, 'device_used_bytes': device_used_bytes,
        'system_ceiling_bytes': system_ceiling_bytes,
        'process_ceiling_bytes': process_ceiling_bytes,
        'estimated_working_set_bytes': estimate,
        'reserved_headroom_bytes': reserve,
        'rematerialization': bool(rematerialization), 'recycles': recycles,
        'recycle_mode': str(recycle_mode), 'model_kind': model_kind,
        'measured_reference': {'residues': 406, 'peak_process_footprint_gib': 27.48,
                               'runtime': 'metaljax 0.11.9+fbc3',
                               'note': 'Two processes, one BF16 multimer fixture on M5; not a memory guarantee'},
        'note': ('The estimate includes a safety margin and is not a guarantee; '
                 'the native runtime guards remain authoritative.'),
    }
    if total_bytes and reserve is not None and score_bytes > max(0, total_bytes - reserve):
        raise MemoryError(f'Unchunked Metal attention needs at least {score_bytes / _GIB:.2f} GiB '
                          'for one score tensor, exceeding physical memory after system headroom; '
                          'reduce the target size or use a larger-memory host.')
    if device_limit_bytes and score_bytes > device_limit_bytes:
        raise MemoryError(f'Unchunked Metal attention needs at least {score_bytes / _GIB:.2f} GiB '
                          f'for one score tensor, exceeding the {device_limit_bytes / _GIB:.2f} GiB '
                          'device limit; reduce the target size or use a larger-memory host.')
    if process_ceiling_bytes and score_bytes > process_ceiling_bytes:
        raise MemoryError('One unchunked Metal score tensor exceeds the configured process limit; '
                          'reduce the target size or use a larger-memory host.')

    # The total process estimate and per-device allocation limit are different
    # metrics. Do not compare the former against device remaining capacity.
    # A free-device comparison can be inconclusive when buffers may be reused.
    uncertain = []
    if device_limit_bytes and device_used_bytes is not None:
        if score_bytes > max(0, device_limit_bytes - device_used_bytes):
            uncertain.append('device buffers may already include reusable storage')
    if estimate is None:
        result['note'] = ('Whole-model peak is uncalibrated for this rematerialization '
                          'or recycle mode; native runtime guards remain authoritative.')
        result['uncertainty'] = uncertain
        return result

    # A known retained allocation overlaps the estimate; all other current
    # footprint remains additive unless separately proven reclaimable. The
    # optimistic bound below is used only to distinguish refusal from unknown,
    # never to claim that unclassified memory will actually be available.
    upper_extra = max(0, estimate - retained)
    possible_overlap = min(upper_extra, unclassified) if unclassified is not None else upper_extra
    constraints = []
    if available_bytes is not None and reserve is not None:
        constraints.append(('system headroom', available_bytes - reserve, 0, upper_extra, possible_overlap))
    if system_ceiling_bytes:
        if system_claimed_bytes is None:
            uncertain.append('system claimed memory is unavailable')
        else:
            system_retained = min(retained, system_claimed_bytes)
            system_extra = max(0, estimate - system_retained)
            system_unknown = (min(unclassified, max(0, system_claimed_bytes - system_retained))
                              if unclassified is not None else system_extra)
            constraints.append(('system claimed-memory limit',
                                system_ceiling_bytes - system_claimed_bytes, claimed_cache,
                                system_extra, min(system_extra, system_unknown)))
    if process_ceiling_bytes:
        if process_bytes is None:
            uncertain.append('process physical footprint is unavailable')
        else:
            constraints.append(('process physical-footprint limit',
                                process_ceiling_bytes - process_bytes, cache, upper_extra, possible_overlap))
    if not constraints:
        uncertain.append('no comparable memory headroom is available')
    headrooms = []
    violations = []
    for label, headroom, confirmed_cache, required, unknown_overlap in constraints:
        extra = required - confirmed_cache
        optimistic = extra - min(required, unknown_overlap)
        headrooms.append(headroom)
        if optimistic > headroom:
            violations.append((label, optimistic, headroom))
        elif extra > headroom:
            uncertain.append(label + ' depends on unclassified live allocations or cache')
    result['estimated_additional_bytes'] = upper_extra
    result['possible_unclassified_overlap_bytes'] = possible_overlap
    result['projected_process_footprint_bytes'] = (
        process_bytes + upper_extra - cache if process_bytes is not None else None)
    result['projected_system_claimed_bytes'] = (
        system_claimed_bytes + max(0, estimate - min(retained, system_claimed_bytes)) - claimed_cache
        if system_claimed_bytes is not None else None)
    if headrooms:
        result['usable_bytes'] = max(0, min(headrooms))
    result['uncertainty'] = uncertain
    if violations:
        label, needed, headroom = violations[0]
        raise MemoryError(f'The estimated unchunked Metal working set is {estimate / _GIB:.2f} GiB; '
                          f'at least {max(0, needed) / _GIB:.2f} GiB of additional headroom is estimated '
                          f'but {label} has {max(0, headroom) / _GIB:.2f} GiB. This estimate includes a '
                          'safety margin and is not a guarantee; reduce the target size or use a '
                          'larger-memory host.')
    result['status'] = 'unknown' if uncertain else 'admitted'
    return result


def _process_physical_footprint():
    """Return the same process physical-footprint metric used by MetalJAX."""
    import ctypes
    import sys
    if sys.platform != 'darwin':
        return None
    class RusageInfoV4(ctypes.Structure):
        _fields_ = [('ri_uuid', ctypes.c_uint8 * 16)] + [
            ('ri_' + name, ctypes.c_uint64) for name in (
                'user_time system_time pkg_idle_wkups interrupt_wkups pageins '
                'wired_size resident_size phys_footprint proc_start_abstime proc_exit_abstime '
                'child_user_time child_system_time child_pkg_idle_wkups child_interrupt_wkups '
                'child_pageins child_elapsed_abstime diskio_bytesread diskio_byteswritten '
                'cpu_time_qos_default cpu_time_qos_maintenance cpu_time_qos_background '
                'cpu_time_qos_utility cpu_time_qos_legacy cpu_time_qos_user_initiated '
                'cpu_time_qos_user_interactive billed_system_time serviced_system_time '
                'logical_writes lifetime_max_phys_footprint instructions cycles billed_energy '
                'serviced_energy interval_max_phys_footprint runnable_time').split()]
    try:
        lib = ctypes.CDLL('/usr/lib/libproc.dylib', use_errno=True)
        lib.proc_pid_rusage.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
        lib.proc_pid_rusage.restype = ctypes.c_int
        usage = RusageInfoV4()
        if lib.proc_pid_rusage(os.getpid(), 4, ctypes.byref(usage)) == 0:
            return int(usage.ri_phys_footprint)
    except (OSError, AttributeError):
        pass
    return None


def _system_claimed_memory():
    """Read MetalJAX's host_statistics64 claimed formula via macOS vm_stat."""
    import re
    import subprocess
    import sys
    if sys.platform != 'darwin':
        return None
    try:
        env = dict(os.environ, LC_ALL='C')
        output = subprocess.run(['/usr/bin/vm_stat'], check=True, capture_output=True,
                                text=True, timeout=3, env=env).stdout
        match = re.search(r'page size of (\d+) bytes', output)
        pages = {m.group(1): int(m.group(2)) for line in output.splitlines()
                 if (m := re.match(r'^([^:]+):\s*(\d+)\.', line))}
        names = ('Pages wired down', 'Anonymous pages', 'Pages purgeable',
                 'Pages occupied by compressor')
        if not match or not all(name in pages for name in names):
            return None
        return int(match.group(1)) * (pages[names[0]] + max(0, pages[names[1]] - pages[names[2]]) + pages[names[3]])
    except (OSError, subprocess.SubprocessError, ValueError):
        return None


def _verified_retained_model_pages(model_params):
    """Count only fully covered, resident, dirty anonymous model buffer pages.

    Called only for an already prepared model. JAX's Metal external reference
    may settle an array, so a first preparation never probes model buffers.
    Overlapping ranges are merged; incomplete pages and any failed queries
    receive no credit. This intentionally does not count a logical array size as
    a resident allocation or claim unrelated process memory as reusable.
    """
    import ctypes
    import sys
    if model_params is None or sys.platform != 'darwin':
        return None, {'source': 'unavailable or first preparation'}
    try:
        import jax
        import numpy as np
        page = int(os.sysconf('SC_PAGE_SIZE'))
        lib = ctypes.CDLL(None, use_errno=True)
        task = ctypes.c_uint.in_dll(lib, 'mach_task_self_').value
        query = lib.mach_vm_page_range_query
        query.argtypes = [ctypes.c_uint, ctypes.c_uint64, ctypes.c_uint64,
                          ctypes.c_uint64, ctypes.POINTER(ctypes.c_uint64)]
        query.restype = ctypes.c_int
    except (ImportError, AttributeError, OSError, ValueError):
        return None, {'source': 'resident-page query unavailable'}
    stack, seen, ranges, failures = [model_params], set(), [], 0
    checked = 0
    while stack:
        value = stack.pop()
        if id(value) in seen:
            continue
        seen.add(id(value))
        if isinstance(value, dict):
            stack.extend(value.values())
            continue
        if isinstance(value, (list, tuple)):
            stack.extend(value)
            continue
        try:
            if isinstance(value, np.ndarray):
                if not value.flags.owndata or not value.flags.c_contiguous or value.dtype.hasobject or isinstance(value, np.memmap):
                    continue
                ptr, size = int(value.ctypes.data), int(value.nbytes)
            elif isinstance(value, jax.Array):
                device = value.device
                if callable(device):
                    device = device()
                if getattr(device, 'platform', None) != 'metal' or not value.is_fully_addressable:
                    continue
                ptr, size = int(value.unsafe_buffer_pointer()), int(value.nbytes)
            else:
                continue
            checked += 1
            start = ((ptr + page - 1) // page) * page
            end = ((ptr + size) // page) * page
            if start < end:
                ranges.append((start, end))
        except (AttributeError, RuntimeError, ValueError, TypeError):
            failures += 1
    merged = []
    for start, end in sorted(ranges):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    good_pages = queried_pages = 0
    # Exclude pages that could belong to an external mapping or cache, or have
    # been paged out; only resident dirty private-like model pages are counted.
    excluded = 0x0002 | 0x0010 | 0x0020 | 0x0040 | 0x0080 | 0x0800
    for start, end in merged:
        while start < end:
            count = min((end - start) // page, 16384)
            flags = (ctypes.c_int * count)()
            returned = ctypes.c_uint64(count)
            if query(task, start, count * page, ctypes.addressof(flags), ctypes.byref(returned)) == 0:
                n = min(count, int(returned.value))
                queried_pages += n
                good_pages += sum((f & 0x0009) == 0x0009 and not (f & excluded) for f in flags[:n])
            else:
                failures += 1
            start += count * page
    return good_pages * page, {
        'source': 'resident dirty anonymous full pages of retained model buffers',
        'buffers_checked': checked, 'pages_queried': queried_pages,
        'pages_credited': good_pages, 'query_failures': failures,
    }


def current_metal_memory(total_residues, num_heads, element_bytes, *,
                         device_limit_bytes=None, device_used_bytes=None,
                         reported_device_cache_bytes=None, model_params=None,
                         rematerialization=True, recycles=None, recycle_mode='last', model_kind=None):
    # Probe retained pages before the memory samples, so any settling of model
    # arrays is included in the footprint and available-memory measurements.
    retained, evidence = _verified_retained_model_pages(model_params)
    try:
        import psutil
        memory = psutil.virtual_memory()
        total, available = memory.total, memory.available
    except (ImportError, OSError, AttributeError):
        total, available = None, None
    process, claimed = _process_physical_footprint(), _system_claimed_memory()

    def ceiling(name):
        value = os.environ.get(name)
        try:
            mb = int(value) if value else 0
        except ValueError:
            return None
        return mb * 1024 ** 2 if mb > 0 else None

    result = assess_metal_memory(
        total_residues, num_heads, element_bytes, total_bytes=total,
        available_bytes=available, process_bytes=process, system_claimed_bytes=claimed,
        retained_model_bytes=retained, reported_device_cache_bytes=reported_device_cache_bytes,
        device_limit_bytes=device_limit_bytes, device_used_bytes=device_used_bytes,
        system_ceiling_bytes=ceiling('METALJAX_MEM_SYS_MB'),
        process_ceiling_bytes=ceiling('METALJAX_MEM_BUDGET_MB'),
        rematerialization=rematerialization, recycles=recycles,
        recycle_mode=recycle_mode, model_kind=model_kind)
    result['retained_model_measurement'] = evidence
    return result
