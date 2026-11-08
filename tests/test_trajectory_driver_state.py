"""Exercise the real driver loop with CPU-only scientific stubs."""
import ast
import contextlib
import copy
import gc
import importlib.util
import io
import json
import os
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
import warnings

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('driver_scheduler_test_module', ROOT/'functions/trajectory_scheduler.py')
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
source = ast.parse((ROOT/'bindcraft.py').read_text())
start = next(i for i, n in enumerate(source.body)
             if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'script_start_time' for t in n.targets))
end = next(i for i, n in enumerate(source.body)
           if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'elapsed_time' for t in n.targets))
DRIVER_LOOP = compile(ast.fix_missing_locations(ast.Module(body=source.body[start:end], type_ignores=[])),
                      str(ROOT/'bindcraft.py'), 'exec')


class FakeRuntime:
    def __init__(self):
        self.model = object()
    def release_outputs(self):
        pass
    def invalidate(self):
        self.model = None


class FakeTrajectory:
    def __init__(self, length, terminated):
        self.length = length
        self.aux = {'log': {'terminate': terminated}}
        self._tmp = {'best': {'aux': {'log': {'plddt': 0.8, 'ptm': 0.7, 'i_ptm': 0.6,
                                             'pae': 0.1, 'i_pae': 0.2}}}}
    def get_seq(self, get_best=True):
        return ['A' * self.length]


class FakeModel:
    def prep_inputs(self, **kwargs):
        pass


class DriverHarness:
    def __init__(self, root):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.pdb = self.root/'target.pdb'
        self.pdb.write_text('TARGET')
        keys = ['Accepted', 'Accepted/Ranked', 'Trajectory', 'Trajectory/Relaxed',
                'Trajectory/LowConfidence', 'Trajectory/Clashing', 'Trajectory/Plots',
                'Trajectory/Animation', 'MPNN', 'MPNN/Binder', 'MPNN/Relaxed', 'Rejected']
        self.paths = {key:str(self.root/key) for key in keys}
        for path in self.paths.values():
            Path(path).mkdir(parents=True, exist_ok=True)
        self.counter = 0
        self.all_lengths = []
        self.prediction_calls = []
        self.advanced = dict(design_algorithm='4stage', enable_mpnn=True, num_seqs=1,
            force_reject_AA=False, omit_AAs='', optimise_beta=False, num_recycles_validation=3,
            af_params_dir='unused', predict_initial_guess=False, predict_bigbang=False,
            rm_template_seq_predict=False, rm_template_sc_predict=False,
            save_mpnn_fasta=False, remove_unrelaxed_complex=False, remove_binder_monomer=False,
            save_design_animations=False, remove_unrelaxed_trajectory=False,
            max_mpnn_sequences=1, start_monitoring=100, enable_rejection_check=False)
        self.target = dict(design_path=str(self.root), binder_name='UBQ', starting_pdb=str(self.pdb),
                           chains='A', target_hotspot_residues='A44', lengths=[65,200],
                           number_of_final_designs=20)

    def run(self, outcomes, *, goal=20, max_attempts=None, stop_after=None, cache=True,
            attempts_per_length=5, accepted_per_length=1):
        outcomes = iter(outcomes)
        self.target['number_of_final_designs'] = goal
        calls = []
        self.finalization_calls = []
        current = {'outcome': None, 'length': None}
        outer = self
        class Random:
            def randint(self, *args, **kwargs):
                outer.counter += 1
                return [outer.counter]
            def choice(self, values):
                return values[len(values)//2]
        class Shutil:
            def copy(self, source, dest):
                (Path(dest)/Path(source).name).write_text('accepted')
        def binder(name, pdb, chains, hotspot, length, seed, helicity, models, settings, paths, failure, runtime=None):
            result = next(outcomes)
            current.update(outcome=result, length=length)
            calls.append((length, result))
            self.all_lengths.append(length)
            if result == 'partial_crash':
                (Path(paths['Accepted'])/(name+'_mpnn1_model1.pdb')).write_text('accepted')
                raise RuntimeError('interrupted')
            if result == 'crash':
                raise RuntimeError('interrupted')
            terminated = {'clash':'Clashing', 'low_confidence':'LowConfidence'}.get(result, '')
            return FakeTrajectory(length, terminated)
        scores = {key:0 for key in ['binder_score','surface_hydrophobicity','interface_sc','interface_packstat',
            'interface_dG','interface_dSASA','interface_dG_SASA_ratio','interface_fraction',
            'interface_hydrophobicity','interface_nres','interface_interface_hbonds',
            'interface_hbond_percentage','interface_delta_unsat_hbonds','interface_delta_unsat_hbonds_percentage']}
        def interface(*args, **kwargs):
            return scores, {}, '' if current['outcome'] == 'no_interface' else 'B1'
        def check_goal(*args, **kwargs):
            if kwargs.get('finalize'):
                self.finalization_calls.append(kwargs['rank_by'])
            return len(list((self.root/'Accepted').glob('*.pdb'))) >= goal
        def check_trajectories(*args, **kwargs):
            return stop_after is not None and len(calls) >= stop_after
        def mpnn(*args, **kwargs):
            count = self.advanced['num_seqs']
            return {'seq':['A' * (current['length']-1) + 'CDEFGH'[i] for i in range(count)],
                    'score':[float(i) for i in range(count)], 'seqid':[0.0] * count}
        def predict_complex(*args, **kwargs):
            self.prediction_calls.append(args[1])
            return ({1:dict(metrics)}, current['outcome'] != 'early_reject', {})
        class FailureFrame:
            columns = []
            def to_csv(self, *args, **kwargs):
                pass
        metrics = {'pLDDT':80.0, 'pTM':0.5, 'i_pTM':0.5, 'pAE':1.0, 'i_pAE':1.0}
        env = dict(time=__import__('time'), os=os, gc=gc, warnings=warnings,
            np=SimpleNamespace(random=Random(), arange=lambda start, end:list(range(start, end))),
            shutil=Shutil(), copy_dict=copy.deepcopy,
            args=SimpleNamespace(rank_by='i_pTM', max_attempts=max_attempts,
                                 attempts_per_length=attempts_per_length,
                                 accepted_per_length=accepted_per_length),
            cache_enabled=cache, TrajectoryDesignRuntime=FakeRuntime,
            TrajectoryScheduleCheckpoint=module.TrajectoryScheduleCheckpoint,
            scheduling_identity=module.scheduling_identity, target_settings=copy.deepcopy(self.target),
            advanced_settings=copy.deepcopy(self.advanced), filters={}, design_paths=self.paths,
            mpnn_csv=str(self.root/'mpnn.csv'), trajectory_csv=str(self.root/'trajectories.csv'),
            final_csv=str(self.root/'final.csv'), failure_csv=str(self.root/'failures.csv'),
            rejected_mpnn_full_stats_csv=str(self.root/'rejected.csv'), filter_column_names_for_rejected_log=[],
            final_labels=[], design_labels=[], settings_file='target', filters_file='filters', advanced_file='advanced',
            design_models=[0], prediction_models=[0], multimer_validation=False, use_pyrosetta=False,
            check_accepted_designs=check_goal, check_n_trajectories=check_trajectories,
            load_helicity=lambda settings:0, binder_hallucination=binder, pr_relax=lambda *a,**k:None,
            calculate_clash_score=lambda *a:0, calc_ss_percentage=lambda *a:(0,0,0,0,0,0,0,0),
            score_interface=interface, validate_design_sequence=lambda *a:'', target_pdb_rmsd=lambda *a:0,
            insert_data=lambda *a:None, mpnn_gen_sequence=mpnn, mk_afdesign_model=lambda **k:FakeModel(),
            predict_binder_complex=predict_complex,
            predict_binder_alone=lambda *a,**k:{1:dict(metrics)},
            calculate_averages=lambda *a,**k:dict(metrics),
            check_filters=lambda *a: {} if current['outcome'] == 'final_reject' else True,
            pd=SimpleNamespace(read_csv=lambda *a,**k:FailureFrame()),
            clear_mem=lambda:None, clear_dssp_cache=lambda:None, rmsd_binder=0, mpnn_interface_residues='B1')
        with contextlib.redirect_stdout(io.StringIO()):
            exec(DRIVER_LOOP, env)
        return calls


class DriverStateTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def test_driver_resumes_next_region_after_completed_success_block(self):
        resumed = DriverHarness(Path(self.tmp.name)/'resumed')
        first = resumed.run(['success'], stop_after=1)
        second = resumed.run(['success'], max_attempts=1)
        uninterrupted = DriverHarness(Path(self.tmp.name)/'uninterrupted')
        whole = uninterrupted.run(['success']*2, stop_after=2)
        self.assertEqual([n for n,_ in first+second], [n for n,_ in whole])
        self.assertEqual([n for n,_ in first+second], [82,116])
        self.assertEqual(len(second), 1)  # --max-attempts counts only this invocation

    def test_driver_flushes_failure_and_no_interface_continue_before_attempt_limit(self):
        for outcome in ('clash', 'low_confidence', 'no_interface'):
            with self.subTest(outcome=outcome):
                h = DriverHarness(Path(self.tmp.name)/outcome)
                calls = h.run([outcome, outcome], goal=5, max_attempts=2)
                self.assertEqual(len(calls), 2)
                state = json.loads((h.root/'.trajectory-schedule.json').read_text())
                self.assertIsNone(state['inflight'])
                self.assertEqual(sum(b['attempted'] for b in state['scheduler']['bins']), 2)
                self.assertEqual(sum(b['accepted'] for b in state['scheduler']['bins']), 0)

    def test_driver_failure_block_resumes_at_the_same_next_region(self):
        resumed = DriverHarness(Path(self.tmp.name)/'failed-resumed')
        first = resumed.run(['clash']*5, stop_after=5)
        second = resumed.run(['clash'], max_attempts=1)
        uninterrupted = DriverHarness(Path(self.tmp.name)/'failed-whole')
        whole = uninterrupted.run(['clash']*6, stop_after=6)
        self.assertEqual([n for n,_ in first+second], [n for n,_ in whole])
        self.assertEqual([n for n,_ in first+second], [82]*5 + [116])

    def test_driver_honors_configurable_repeats_and_still_stops_on_success(self):
        h = DriverHarness(Path(self.tmp.name)/'custom-repeats')
        calls = h.run(['clash']*2 + ['success', 'clash'], max_attempts=4, attempts_per_length=2)
        self.assertEqual([n for n, _ in calls], [82, 82, 116, 150])
        state = json.loads((h.root/'.trajectory-schedule.json').read_text())
        self.assertEqual(state['attempts_per_visit'], 2)

        h = DriverHarness(Path(self.tmp.name)/'long-repeats')
        calls = h.run(['clash']*10, max_attempts=10, attempts_per_length=9)
        self.assertEqual([n for n, _ in calls], [82]*9 + [116])

    def test_driver_continues_until_requested_accepted_designs(self):
        h = DriverHarness(Path(self.tmp.name)/'accepted-target')
        calls = h.run(['success', 'clash', 'success', 'success', 'clash'],
                      max_attempts=5, attempts_per_length=5, accepted_per_length=3)
        self.assertEqual([n for n, _ in calls], [82, 82, 82, 82, 116])
        state = json.loads((h.root/'.trajectory-schedule.json').read_text())
        self.assertEqual(state['accepted_per_visit'], 3)
        self.assertEqual(state['scheduler']['bins'][0]['accepted'], 3)

    def test_driver_attempt_limit_can_end_visit_before_acceptance_target(self):
        h = DriverHarness(Path(self.tmp.name)/'accepted-target-attempt-limit')
        calls = h.run(['success', 'clash', 'clash'], max_attempts=3,
                      attempts_per_length=2, accepted_per_length=2)
        self.assertEqual([n for n, _ in calls], [82, 82, 116])

    def test_driver_counts_multiple_accepted_designs_in_one_trajectory(self):
        h = DriverHarness(Path(self.tmp.name)/'accepted-target-multiple')
        h.advanced.update(num_seqs=3, max_mpnn_sequences=3)
        calls = h.run(['success', 'clash'], goal=10, max_attempts=2,
                      accepted_per_length=2)
        self.assertEqual([n for n, _ in calls], [82, 116])
        self.assertEqual(len(list((h.root/'Accepted').glob('*.pdb'))), 3)
        self.assertEqual(len(h.prediction_calls), 3)

    def test_driver_resumes_with_a_higher_accepted_design_target(self):
        h = DriverHarness(Path(self.tmp.name)/'accepted-target-resume')
        self.assertEqual([n for n, _ in h.run(['success'], max_attempts=1)], [82])
        with self.assertWarnsRegex(RuntimeWarning, 'accepted-design target changed from 1 to 3'):
            calls = h.run(['success', 'success', 'clash'], max_attempts=3,
                          accepted_per_length=3)
        self.assertEqual([n for n, _ in calls], [82, 82, 116])
        state = json.loads((h.root/'.trajectory-schedule.json').read_text())
        self.assertEqual(state['scheduler']['bins'][0]['accepted'], 3)

    def test_driver_applies_changed_limit_when_resuming(self):
        h = DriverHarness(Path(self.tmp.name)/'changed-repeats')
        self.assertEqual([n for n, _ in h.run(['clash']*3, max_attempts=3)], [82]*3)
        with self.assertWarnsRegex(RuntimeWarning, 'repeat limit changed from 5 to 2'):
            calls = h.run(['clash'], max_attempts=1, attempts_per_length=2)
        self.assertEqual([n for n, _ in calls], [116])

    def test_driver_mid_block_resume_matches_uninterrupted_progress(self):
        resumed = DriverHarness(Path(self.tmp.name)/'middle-resumed')
        first = resumed.run(['clash']*3, stop_after=3)
        second = resumed.run(['clash','success','success'], stop_after=3)
        uninterrupted = DriverHarness(Path(self.tmp.name)/'middle-whole')
        whole = uninterrupted.run(['clash']*4+['success','success'], stop_after=6)
        self.assertEqual([n for n,_ in first+second], [n for n,_ in whole])
        self.assertEqual([n for n,_ in first+second], [82]*5+[116])

    def test_preliminary_success_is_not_final_filter_acceptance(self):
        for rejection in ('early_reject', 'final_reject'):
            with self.subTest(rejection=rejection):
                h = DriverHarness(Path(self.tmp.name)/rejection)
                calls = h.run([rejection, 'success', 'clash'], max_attempts=3)
                self.assertEqual([n for n,_ in calls], [82,82,116])
                state = json.loads((h.root/'.trajectory-schedule.json').read_text())
                self.assertEqual(state['scheduler']['bins'][0]['attempted'], 2)
                self.assertEqual(state['scheduler']['bins'][0]['accepted'], 1)

    def test_success_on_first_middle_and_fifth_attempt_moves_to_next_bin(self):
        for success_at in (1, 3, 5):
            with self.subTest(success_at=success_at):
                h = DriverHarness(Path(self.tmp.name)/f'success-{success_at}')
                calls = h.run(['clash']*(success_at-1)+['success','clash'],
                              max_attempts=success_at+1)
                self.assertEqual([n for n,_ in calls], [82]*success_at+[116])

    def test_multiple_accepted_outputs_respect_mpnn_limit_and_global_goal(self):
        for cache in (True, False):
            with self.subTest(cache=cache):
                h = DriverHarness(Path(self.tmp.name)/f'multiple-{cache}')
                h.advanced.update(num_seqs=3, max_mpnn_sequences=3)
                calls = h.run(['success'], goal=2, cache=cache)
                self.assertEqual(len(calls), 1)
                self.assertEqual(len(list((h.root/'Accepted').glob('*.pdb'))), 2)
                self.assertEqual(len(h.prediction_calls), 2)
        h = DriverHarness(Path(self.tmp.name)/'mpnn-cap')
        h.advanced.update(num_seqs=3, max_mpnn_sequences=2)
        calls = h.run(['success','clash'], goal=10, max_attempts=2)
        self.assertEqual([n for n,_ in calls], [82,116])
        self.assertEqual(len(list((h.root/'Accepted').glob('*.pdb'))), 2)
        self.assertEqual(len(h.prediction_calls), 2)

    def test_driver_changed_goal_and_range_discard_old_history(self):
        h = DriverHarness(Path(self.tmp.name)/'changed')
        h.run(['success'], goal=20, stop_after=1)
        with self.assertWarnsRegex(RuntimeWarning, 'campaign settings changed'):
            h.run(['success'], goal=21, stop_after=1)
        state = json.loads((h.root/'.trajectory-schedule.json').read_text())
        self.assertEqual(sum(b['attempted'] for b in state['scheduler']['bins']), 1)
        h.target['lengths'] = [100, 103]
        with self.assertWarnsRegex(RuntimeWarning, 'campaign settings changed'):
            self.assertEqual(h.run(['clash'], goal=21, stop_after=1)[0][0], 100)
        state = json.loads((h.root/'.trajectory-schedule.json').read_text())
        self.assertEqual(sum(b['attempted'] for b in state['scheduler']['bins']), 1)

    def test_cache_off_does_not_create_or_update_checkpoint(self):
        h = DriverHarness(Path(self.tmp.name)/'off')
        h.run(['clash'], goal=4, max_attempts=1, cache=False)
        self.assertFalse((h.root/'.trajectory-schedule.json').exists())
        h.run(['clash'], goal=4, max_attempts=1)
        before = (h.root/'.trajectory-schedule.json').read_bytes()
        h.run(['clash'], goal=4, max_attempts=1, cache=False)
        self.assertEqual((h.root/'.trajectory-schedule.json').read_bytes(), before)

    def test_attempt_cap_finalizes_partial_and_empty_results_with_cache_on_or_off(self):
        for cache in (True, False):
            for outcome in ('success', 'clash'):
                with self.subTest(cache=cache, outcome=outcome):
                    h = DriverHarness(Path(self.tmp.name)/f'partial-{cache}-{outcome}')
                    calls = h.run([outcome], goal=10, max_attempts=1, cache=cache)
                    self.assertEqual(len(calls), 1)
                    self.assertEqual(h.finalization_calls, ['Average_i_pTM'])

    def test_driver_completion_flushes_last_accepted_attempt(self):
        h = DriverHarness(Path(self.tmp.name)/'completion')
        calls = h.run(['success','success'], goal=2)
        self.assertEqual(len(calls), 2)
        self.assertEqual(h.finalization_calls, [])
        state = json.loads((h.root/'.trajectory-schedule.json').read_text())
        self.assertIsNone(state['inflight'])
        self.assertEqual(len(state['accepted_files']), 2)
        self.assertEqual(sum(b['successful'] for b in state['scheduler']['bins']), 2)

    def test_driver_recovers_interruption_without_spending_new_attempt_budget(self):
        h = DriverHarness(Path(self.tmp.name)/'crash')
        with self.assertRaisesRegex(RuntimeError, 'interrupted'):
            h.run(['crash'], goal=4, max_attempts=4)
        state = json.loads((h.root/'.trajectory-schedule.json').read_text())
        self.assertIsNotNone(state['inflight'])
        with self.assertWarnsRegex(RuntimeWarning, 'Recovered an interrupted'):
            calls = h.run(['clash'], goal=4, max_attempts=1)
        self.assertEqual(len(calls), 1)
        state = json.loads((h.root/'.trajectory-schedule.json').read_text())
        self.assertIsNone(state['inflight'])
        self.assertEqual(sum(b['attempted'] for b in state['scheduler']['bins']), 2)

    def test_driver_recovers_accepted_file_from_interrupted_attempt(self):
        h = DriverHarness(Path(self.tmp.name)/'partial-crash')
        with self.assertRaisesRegex(RuntimeError, 'interrupted'):
            h.run(['partial_crash'], goal=4, max_attempts=4)
        with self.assertWarnsRegex(RuntimeWarning, 'Recovered an interrupted'):
            calls = h.run(['clash'], goal=4, max_attempts=1)
        self.assertEqual([length for length, _ in calls], [116])
        state = json.loads((h.root/'.trajectory-schedule.json').read_text())
        self.assertEqual(len(state['accepted_files']), 1)
        self.assertEqual(sum(b['successful'] for b in state['scheduler']['bins']), 1)
        self.assertEqual(sum(b['attempted'] for b in state['scheduler']['bins']), 2)

    def test_driver_starts_legacy_output_from_existing_files(self):
        h = DriverHarness(Path(self.tmp.name)/'legacy')
        (h.root/'Accepted'/'old_mpnn1_model1.pdb').write_text('accepted')
        calls = h.run(['clash'], goal=20, max_attempts=1)
        self.assertEqual(len(calls), 1)
        state = json.loads((h.root/'.trajectory-schedule.json').read_text())
        self.assertEqual(state['accepted_files'], ['old_mpnn1_model1.pdb'])
        self.assertEqual(sum(b['attempted'] for b in state['scheduler']['bins']), 1)


if __name__ == '__main__':
    unittest.main()
