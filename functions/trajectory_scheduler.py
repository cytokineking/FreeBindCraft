"""Exact-length scheduling and crash-safe checkpoints for the optional cache.

The scheduler is independent of JAX. Its checkpoint reconciles progress with
the accepted filenames reported by the driver.
"""
from dataclasses import dataclass, field
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile
import warnings

DEFAULT_TRAJECTORY_REPEATS = 5
DEFAULT_TRAJECTORY_ACCEPTED_DESIGNS = 1


@dataclass
class _Bin:
    low: int
    high: int
    attempted: int = 0
    successful: int = 0
    accepted: int = 0
    visited: set = field(default_factory=set)

    def next_length(self):
        # First cover the middle, then repeatedly choose the integer furthest
        # from the already sampled lengths. Ties go to the smaller length.
        if len(self.visited) == self.high - self.low + 1:
            self.visited.clear()
        if not self.visited:
            length = (self.low + self.high + 1) // 2
        else:
            length = max((n for n in range(self.low, self.high + 1)
                          if n not in self.visited),
                         key=lambda n: (min(abs(n - v) for v in self.visited), -n))
        self.visited.add(length)
        return length


class TrajectoryLengthScheduler:
    """Interleave regions until a visit reaches its accepted-design or attempt limit."""

    def __init__(self, min_length, max_length, remaining_designs, max_attempts=None,
                 *, attempts_per_visit=DEFAULT_TRAJECTORY_REPEATS,
                 accepted_per_visit=DEFAULT_TRAJECTORY_ACCEPTED_DESIGNS):
        if min_length > max_length or remaining_designs < 1:
            raise ValueError("invalid length range or remaining design count")
        if max_attempts is not None and max_attempts < 1:
            raise ValueError("max_attempts must be positive")
        count = min(4, max_length - min_length + 1, remaining_designs)
        width, extra = divmod(max_length - min_length + 1, count)
        self.bins = []
        low = min_length
        for index in range(count):
            size = width + (index < extra)
            self.bins.append(_Bin(low, low + size - 1))
            low += size
        self._next_bin = 0
        self.current = None
        self.set_visit_limits(attempts_per_visit, accepted_per_visit)

    def set_visit_limits(self, attempts_per_visit, accepted_per_visit):
        """Apply new limits, including to a previously started visit."""
        if type(attempts_per_visit) is not int or attempts_per_visit < 1:
            raise ValueError('attempts_per_visit must be a positive integer')
        if type(accepted_per_visit) is not int or accepted_per_visit < 1:
            raise ValueError('accepted_per_visit must be a positive integer')
        self.attempts_per_visit = attempts_per_visit
        self.accepted_per_visit = accepted_per_visit
        if self.current is not None:
            self.current['done'] = (self.current['accepted'] >= accepted_per_visit or
                                    self.current['attempts'] >= attempts_per_visit)

    def next_length(self, remaining_designs, remaining_attempts=None):
        if remaining_designs < 1 or (remaining_attempts is not None and remaining_attempts < 1):
            raise ValueError("a visit needs at least one remaining design and attempt")
        if self.current is not None and not self.current['done']:
            return self.current['length']

        index = self._next_bin
        region = self.bins[index]
        self.current = dict(bin=index, length=region.next_length(), attempts=0,
                            successful=0, accepted=0, done=False)
        self._next_bin = (index + 1) % len(self.bins)
        return self.current['length']

    def record(self, accepted_outputs, *, recovering_legacy=False):
        """Record one attempted trajectory, including zero-output failures."""
        if (self.current is None or (self.current['done'] and not recovering_legacy) or
                type(accepted_outputs) is not int or accepted_outputs < 0):
            raise ValueError("no active visit or invalid accepted count")
        block = self.current
        region = self.bins[block['bin']]
        block['attempts'] += 1
        region.attempted += 1
        block['accepted'] += accepted_outputs
        region.accepted += accepted_outputs
        if accepted_outputs:
            block['successful'] += 1
            region.successful += 1
        block['done'] = (block['accepted'] >= self.accepted_per_visit or
                         block['attempts'] >= self.attempts_per_visit)

    def to_state(self):
        return {
            'bins': [dict(low=b.low, high=b.high, attempted=b.attempted,
                          successful=b.successful, accepted=b.accepted,
                          visited=sorted(b.visited)) for b in self.bins],
            'next_bin': self._next_bin,
            'current': None if self.current is None else dict(self.current),
        }

    @classmethod
    def from_state(cls, state, min_length, max_length, *, legacy=False,
                   fixed_v2=False, attempts_per_visit=DEFAULT_TRAJECTORY_REPEATS,
                   accepted_per_visit=DEFAULT_TRAJECTORY_ACCEPTED_DESIGNS):
        """Validate saved state before it can affect scheduling."""
        def integer(value, minimum=0):
            if type(value) is not int or value < minimum:
                raise ValueError('invalid scheduler counter')
            return value

        integer(attempts_per_visit, 1)
        integer(accepted_per_visit, 1)
        if not isinstance(state, dict) or not isinstance(state.get('bins'), list):
            raise ValueError('invalid scheduler state')
        data = state['bins']
        if not 1 <= len(data) <= min(4, max_length - min_length + 1):
            raise ValueError('invalid region count')
        width, extra = divmod(max_length - min_length + 1, len(data))
        bins, low = [], min_length
        for index, entry in enumerate(data):
            if not isinstance(entry, dict):
                raise ValueError('invalid region')
            high = low + width + (index < extra) - 1
            if (type(entry.get('low')) is not int or type(entry.get('high')) is not int or
                    entry['low'] != low or entry['high'] != high):
                raise ValueError('region boundaries do not match the campaign')
            attempted = integer(entry.get('attempted'))
            successful = integer(entry.get('successful'))
            accepted = integer(entry.get('accepted'))
            visited = entry.get('visited')
            if (not isinstance(visited, list) or
                    any(type(n) is not int or n < low or n > high for n in visited) or
                    len(visited) != len(set(visited)) or
                    successful > attempted or successful > accepted or
                    bool(successful) != bool(accepted)):
                raise ValueError('invalid region history')
            bins.append(_Bin(low, high, attempted, successful, accepted, set(visited)))
            low = high + 1

        next_bin = integer(state.get('next_bin'))
        if next_bin >= len(bins):
            raise ValueError('invalid next region')
        current = state.get('current')
        if current is None:
            if next_bin != 0 or any(b.attempted or b.visited for b in bins):
                raise ValueError('missing current visit')
        else:
            if not isinstance(current, dict):
                raise ValueError('invalid visit')
            keys = {'bin', 'length', 'attempts', 'successful', 'accepted', 'done'}
            if legacy:
                keys |= {'target', 'limit'}
            if set(current) != keys or type(current['done']) is not bool:
                raise ValueError('invalid visit fields')
            index = integer(current['bin'])
            length = integer(current['length'], 1)
            attempts = integer(current['attempts'])
            successful = integer(current['successful'])
            accepted = integer(current['accepted'])
            if legacy:
                target = integer(current['target'], 1)
                limit = integer(current['limit'], 1)
                if limit > 8 or attempts > limit:
                    raise ValueError('invalid legacy attempt limit')
                expected_done = accepted >= target or attempts >= limit
            else:
                # Version-1 visits may have used six to eight attempts before
                # migration to version 2. Later versions can use any positive limit.
                if fixed_v2 and attempts > 8:
                    raise ValueError('invalid visit attempt count')
                expected_done = accepted >= accepted_per_visit or attempts >= attempts_per_visit
            if index >= len(bins):
                raise ValueError('invalid visit region')
            region = bins[index]
            if (length not in region.visited or
                    successful > attempts or successful > accepted or
                    bool(successful) != bool(accepted) or
                    attempts > region.attempted or successful > region.successful or
                    accepted > region.accepted or next_bin != (index + 1) % len(bins) or
                    current['done'] != expected_done):
                raise ValueError('inconsistent visit history')
            current = {key: current[key] for key in ('bin', 'length', 'attempts', 'successful', 'accepted')}
            current['done'] = accepted >= accepted_per_visit or attempts >= attempts_per_visit
        result = cls.__new__(cls)
        result.bins, result._next_bin, result.current = bins, next_bin, current
        result.attempts_per_visit = attempts_per_visit
        result.accepted_per_visit = accepted_per_visit
        return result


def scheduling_identity(target_settings, advanced_settings, filters):
    """Fingerprint settings and target contents that define a campaign."""
    target = dict(target_settings)
    target.pop('design_path', None)
    pdb = target.get('starting_pdb')
    if isinstance(pdb, str):
        try:
            with open(os.path.expanduser(pdb), 'rb') as stream:
                target['starting_pdb_sha256'] = hashlib.sha256(stream.read()).hexdigest()
        except OSError:
            # The normal input checks still report an unusable target path.
            target['starting_pdb_sha256'] = None
    payload = json.dumps([target, advanced_settings, filters], sort_keys=True,
                         separators=(',', ':'), ensure_ascii=True)
    return hashlib.sha256(payload.encode('utf-8')).hexdigest()


class TrajectoryScheduleCheckpoint:
    """Atomic scheduling history; accepted PDBs remain the source of truth."""

    VERSION = 4

    def __init__(self, path, identity, min_length, max_length,
                 remaining_designs, remaining_attempts, accepted_files,
                 attempts_per_visit=DEFAULT_TRAJECTORY_REPEATS,
                 accepted_per_visit=DEFAULT_TRAJECTORY_ACCEPTED_DESIGNS):
        if type(attempts_per_visit) is not int or attempts_per_visit < 1:
            raise ValueError('attempts_per_visit must be a positive integer')
        if type(accepted_per_visit) is not int or accepted_per_visit < 1:
            raise ValueError('accepted_per_visit must be a positive integer')
        self.path = Path(path)
        self.identity = identity
        self.accepted = self._accepted_names(accepted_files)
        self.inflight = None
        self.scheduler = None
        try:
            with self.path.open(encoding='utf-8') as stream:
                state = json.load(stream)
        except FileNotFoundError:
            state = None
        except (ValueError, UnicodeError):
            state = None
            warnings.warn('Invalid trajectory scheduling state; starting fresh history.', RuntimeWarning)
        if state is not None:
            try:
                if not isinstance(state, dict) or state.get('version') not in (1, 2, 3, self.VERSION):
                    raise ValueError('unknown checkpoint version')
                legacy = state['version'] == 1
                saved_limit = (state.get('attempts_per_visit') if state['version'] >= 3
                               else DEFAULT_TRAJECTORY_REPEATS)
                saved_target = (state.get('accepted_per_visit') if state['version'] == self.VERSION
                                else DEFAULT_TRAJECTORY_ACCEPTED_DESIGNS)
                if type(saved_limit) is not int or saved_limit < 1:
                    raise ValueError('invalid saved trajectory repeat limit')
                if type(saved_target) is not int or saved_target < 1:
                    raise ValueError('invalid saved accepted-design target')
                if state.get('identity') != identity:
                    raise ValueError('campaign settings changed')
                saved_accepted = self._accepted_names(state['accepted_files'])
                scheduler = TrajectoryLengthScheduler.from_state(state['scheduler'], min_length, max_length,
                    legacy=legacy, fixed_v2=state['version'] == 2, attempts_per_visit=saved_limit,
                    accepted_per_visit=saved_target)
                if sum(b.accepted for b in scheduler.bins) > len(saved_accepted):
                    raise ValueError('accepted history exceeds saved files')
                inflight = state.get('inflight')
                if inflight is None:
                    if saved_accepted != self.accepted:
                        raise ValueError('accepted files changed outside the checkpoint')
                else:
                    saved_current = state['scheduler'].get('current')
                    if (scheduler.current is None or saved_current['done'] or
                            not self._valid_attempt_name(inflight, scheduler.current['length'])):
                        raise ValueError('invalid in-flight attempt')
                    # A crash can happen before or after accepted files are copied.
                    # Count only new files attributable to this one trajectory.
                    added = self._added_for_attempt(saved_accepted, self.accepted, inflight)
                    scheduler.record(len(added), recovering_legacy=legacy)
                changed_limit = saved_limit != attempts_per_visit
                changed_target = saved_target != accepted_per_visit
                scheduler.set_visit_limits(attempts_per_visit, accepted_per_visit)
                self.scheduler = scheduler
                if inflight is not None or state['version'] != self.VERSION or changed_limit or changed_target:
                    self._save()  # Commit recovery/migration before a new attempt.
                if inflight is not None:
                    warnings.warn('Recovered an interrupted trajectory as a completed attempt.', RuntimeWarning)
                if legacy:
                    label = 'five' if attempts_per_visit == 5 else str(attempts_per_visit)
                    warnings.warn(f'Migrated trajectory scheduling state to the {label}-attempt policy.', RuntimeWarning)
                else:
                    if changed_limit:
                        warnings.warn(f'Trajectory repeat limit changed from {saved_limit} to {attempts_per_visit}; '
                                      'the current visit uses the new limit.', RuntimeWarning)
                    if changed_target:
                        warnings.warn(f'Trajectory accepted-design target changed from {saved_target} to '
                                      f'{accepted_per_visit}; the current visit uses the new target.', RuntimeWarning)
            except (ValueError, TypeError, KeyError) as error:
                warnings.warn(f'Cannot reuse trajectory scheduling state ({error}); starting fresh history.', RuntimeWarning)
        if self.scheduler is None:
            self.scheduler = TrajectoryLengthScheduler(min_length, max_length,
                remaining_designs, remaining_attempts, attempts_per_visit=attempts_per_visit,
                accepted_per_visit=accepted_per_visit)

    @staticmethod
    def _accepted_names(files):
        if not isinstance(files, (list, tuple, set, frozenset)):
            raise ValueError('invalid accepted-file snapshot')
        if any(not isinstance(name, str) or Path(name).name != name or
               not name.endswith('.pdb') or name.startswith('.') for name in files):
            raise ValueError('invalid accepted-file name')
        if len(files) != len(set(files)):
            raise ValueError('duplicate accepted-file name')
        return set(files)

    @staticmethod
    def _added_for_attempt(before, after, name):
        added = after - before
        pattern = re.compile(re.escape(name) + r'_mpnn\d+_model\d+\.pdb')
        if before - after or any(pattern.fullmatch(item) is None for item in added):
            raise ValueError('accepted files changed outside the current attempt')
        return added

    @staticmethod
    def _valid_attempt_name(name, length):
        if not isinstance(name, str) or Path(name).name != name:
            return False
        match = re.fullmatch(r'.*_l(\d+)_s(\d+)', name)
        return match is not None and int(match.group(1)) == length

    def _save(self):
        state = {'version': self.VERSION, 'identity': self.identity,
                 'attempts_per_visit': self.scheduler.attempts_per_visit,
                 'accepted_per_visit': self.scheduler.accepted_per_visit,
                 'accepted_files': sorted(self.accepted), 'inflight': self.inflight,
                 'scheduler': self.scheduler.to_state()}
        # The temporary file lives beside the checkpoint, so replace is atomic.
        name = None
        try:
            with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=self.path.parent,
                                             prefix='.trajectory-schedule-', delete=False) as stream:
                name = stream.name
                json.dump(state, stream, sort_keys=True, separators=(',', ':'))
                stream.write('\n')
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(name, self.path)
        finally:
            if name is not None and os.path.exists(name):
                os.unlink(name)

    def begin_attempt(self, name, accepted_files):
        """Persist the selected visit before running a trajectory."""
        if self.inflight is not None or self.scheduler.current is None or self.scheduler.current['done']:
            raise ValueError('no active visit for a new attempt')
        if not self._valid_attempt_name(name, self.scheduler.current['length']):
            raise ValueError('invalid trajectory name')
        if self._accepted_names(accepted_files) != self.accepted:
            return False
        self.inflight = name
        self._save()
        return True

    def finish_attempt(self, accepted_files):
        """Commit one attempt, or discard history if external files changed."""
        if self.inflight is None:
            raise ValueError('no in-flight attempt')
        actual = self._accepted_names(accepted_files)
        try:
            added = self._added_for_attempt(self.accepted, actual, self.inflight)
        except ValueError as error:
            warnings.warn(f'{error}; starting fresh trajectory scheduling history.', RuntimeWarning)
            return False
        self.scheduler.record(len(added))
        self.accepted, self.inflight = actual, None
        self._save()
        return True
