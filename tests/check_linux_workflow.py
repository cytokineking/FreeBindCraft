#!/usr/bin/env python3
"""GPU integration checks for installation upgrades and persisted job outputs.

Run in the installed environment:
  python tests/check_linux_workflow.py --case normal --output-dir runs/linux-checks
  python tests/check_linux_workflow.py --case accepted --output-dir runs/linux-checks
  python tests/check_linux_workflow.py --case ranking --output-dir runs/linux-checks

normal uses the unmodified CLI, default filters, and two bounded GPU runs.
accepted replays only the expensive trajectory-generation boundary from a
published PDL1 structure. MPNN, AF2 validation, OpenMM, scoring, CSV persistence
and selection run for real through bindcraft.py with explicit no_filters.
ranking resumes historical accepted outputs through the unmodified CLI,
including CSV migration and an alternate-metric fixture with known ordering.
These checks establish workflow behavior, not new scientifically accepted hits.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'performance_data/pdl1_miniprotein'
REFERENCE_DESIGNS = ['PDL1_l65_s321660_mpnn17', 'PDL1_l65_s321660_mpnn3']

def read_csv(path):
    with path.open(newline='') as stream:
        reader = csv.DictReader(stream)
        return list(reader.fieldnames), list(reader)

def write_csv(path, columns, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def historical_rows():
    columns, rows = read_csv(DATA / 'pdl1_final_design_stats_freebindcraft.csv')
    by_name = {row['Design']: row for row in rows}
    return columns, [by_name[name] for name in REFERENCE_DESIGNS]

def save_reference(design, destination):
    with zipfile.ZipFile(DATA / 'PDL1-freebindcraft.zip') as archive:
        matches = [name for name in archive.namelist()
                   if '/' not in name and re.fullmatch(r'\d+_' + re.escape(design) + r'_model\d\.pdb', name)]
        assert len(matches) == 1, matches
        destination.write_bytes(archive.read(matches[0]))
        return matches[0].split('_', 1)[1]

def configurations(folder, *, target='PDL1', num_seqs=2):
    advanced = json.loads((ROOT / 'settings_advanced/default_4stage_multimer.json').read_text())
    advanced.update(num_seqs=num_seqs, max_mpnn_sequences=2,
                    optimise_beta=False, max_trajectories=False,
                    save_design_trajectory_plots=False, save_design_animations=False,
                    zip_animations=False, zip_plots=False,
                    remove_unrelaxed_trajectory=False, remove_unrelaxed_complex=False,
                    remove_binder_monomer=False)
    settings = {'design_path': str(folder / 'output'), 'binder_name': 'workflow_check',
        'starting_pdb': str(ROOT / ('example/PDL1.pdb' if target == 'PDL1' else 'example/ubiquitin_1ubq.pdb')),
        'chains': 'A', 'target_hotspot_residues': '',
        'lengths': [65, 65] if target == 'PDL1' else [70, 70],
        'number_of_final_designs': 2}
    (folder / 'settings.json').write_text(json.dumps(settings, indent=2))
    (folder / 'advanced.json').write_text(json.dumps(advanced, indent=2))
    return settings

def cli(folder, label, *, rank_by='i_pTM', replay=False, filters='default_filters'):
    if replay:
        command = [sys.executable, str(Path(__file__).resolve()), '--replay-main', str(folder)]
    else:
        command = [sys.executable, str(ROOT / 'bindcraft.py')]
    command += ['--settings', str(folder / 'settings.json'), '--advanced', str(folder / 'advanced.json'),
        '--filters', str(ROOT / f'settings_filters/{filters}.json'), '--no-pyrosetta',
        '--no-plots', '--no-animations', '--verbose', '--seed', '42', '--max-attempts', '1',
        '--rank-by', rank_by]
    log = folder / (label + '.log')
    with log.open('w') as stream:
        result = subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT, timeout=1800)
    assert result.returncode == 0, f'CLI failed; see {log}'
    return log.read_text()

def assert_ranked(folder, expected):
    output = folder / 'output'
    _, final = read_csv(output / 'final_design_stats.csv')
    assert [row['Design'] for row in final] == expected, final
    assert all(row['Rank'].strip() for row in final), 'Accepted designs were not ranked before CLI exit'
    assert [int(float(row['Rank'])) for row in final] == list(range(1, len(expected) + 1)), final
    ranked = sorted((output / 'Accepted/Ranked').glob('*.pdb'))
    assert len(ranked) == len(expected), ranked
    for rank, design in enumerate(expected, 1):
        selected = list((output / 'Accepted/Ranked').glob(f'{rank}_{design}_model*.pdb'))
        accepted = list((output / 'Accepted').glob(f'{design}_model*.pdb'))
        assert len(selected) == len(accepted) == 1
        assert digest(selected[0]) == digest(accepted[0])

def normal(folder):
    configurations(folder, target='ubiquitin', num_seqs=1)
    first_log = cli(folder, 'initial')
    first_names = re.findall(r'^Starting trajectory: (.+)$', first_log, re.M)
    assert len(first_names) == 1, first_names
    output = folder / 'output/Trajectory'
    first_files = {str(p.relative_to(output)): digest(p) for p in output.rglob('*.pdb')}
    assert first_files
    second_log = cli(folder, 'resumed')
    second_names = re.findall(r'^Starting trajectory: (.+)$', second_log, re.M)
    assert len(second_names) == 1 and second_names[0] != first_names[0], second_names
    for relative, before in first_files.items():
        assert digest(output / relative) == before, relative
    assert len({p.name for p in output.rglob('*.pdb')}) == 2
    return {'case': 'normal', 'default_filters': True, 'first_trajectory': first_names[0],
            'resumed_trajectory': second_names[0], 'previous_structures_unchanged': True}

def accepted(folder):
    configurations(folder)
    save_reference(REFERENCE_DESIGNS[0], folder / 'trajectory-reference.pdb')
    log = cli(folder, 'accepted', replay=True, filters='no_filters')
    assert len(re.findall(r'^Starting trajectory: (.+)$', log, re.M)) == 1
    _, rows = read_csv(folder / 'output/mpnn_design_stats.csv')
    assert len(rows) == 2 and len({row['Sequence'] for row in rows}) == 2, rows
    assert log.count('passed all filters') == 2
    expected = [row['Design'] for row in sorted(rows, key=lambda row: float(row['Average_i_pTM']), reverse=True)]
    assert_ranked(folder, expected)
    assert len(list((folder / 'output/MPNN/Binder').glob('*.pdb'))) == 4
    assert len(list((folder / 'output/MPNN/Relaxed').glob('*.pdb'))) == 4
    return {'case': 'accepted', 'trajectory_replayed_from': REFERENCE_DESIGNS[0],
            'real_mpnn_and_af2_validation': True, 'explicit_no_filters': True,
            'workflow_selected_designs': expected, 'scientific_acceptance_claimed': False}

def ranking(folder):
    configurations(folder)
    output = folder / 'output'
    (output / 'Accepted').mkdir(parents=True)
    columns, rows = historical_rows()
    for row in rows:
        temporary = output / 'fixture.pdb'
        name = save_reference(row['Design'], temporary)
        shutil.move(temporary, output / 'Accepted' / name)
    write_csv(output / 'mpnn_design_stats.csv', [c for c in columns if c != 'Rank'], rows)
    write_csv(output / 'final_design_stats.csv', columns, rows)
    preserved = {p.name: digest(p) for p in (output / 'Accepted').glob('*.pdb')}
    initial = cli(folder, 'legacy-resume')
    assert 'Starting trajectory:' not in initial
    assert_ranked(folder, REFERENCE_DESIGNS)
    columns, migrated = read_csv(output / 'mpnn_design_stats.csv')
    assert len(migrated) == 2 and 'Average_ipSAE' in columns
    # Controlled persisted-state values distinguish the two ranking choices.
    scores = dict(zip(REFERENCE_DESIGNS, ['0.20', '0.80']))
    for row in migrated:
        row['Average_ipSAE'] = scores[row['Design']]
    write_csv(output / 'mpnn_design_stats.csv', columns, migrated)
    before_csv = digest(output / 'mpnn_design_stats.csv')
    expected = list(reversed(REFERENCE_DESIGNS))
    for label in ['alternate-ranking', 'repeat-resume']:
        log = cli(folder, label, rank_by='ipSAE')
        assert 'Starting trajectory:' not in log
        assert_ranked(folder, expected)
        assert digest(output / 'mpnn_design_stats.csv') == before_csv
        assert {p.name: digest(p) for p in (output / 'Accepted').glob('*.pdb')} == preserved
    return {'case': 'ranking', 'historical_csv_migrated': True,
            'default_order': REFERENCE_DESIGNS, 'ipSAE_fixture_order': expected,
            'repeat_resume_preserved_inputs': True, 'new_trajectories': 0}

def replay_main(folder, forwarded):
    # The only replaced boundary is expensive de novo trajectory generation.
    # All code after that boundary is the real CLI and real scientific runtime.
    import runpy
    from types import SimpleNamespace
    sys.path.insert(0, str(ROOT))
    import functions
    _, rows = historical_rows()
    row = rows[0]
    metrics = {name: float(row['Average_' + field]) for name, field in
               [('plddt','pLDDT'), ('ptm','pTM'), ('i_ptm','i_pTM'), ('pae','pAE'), ('i_pae','i_pAE')]}
    def replay(design_name, starting_pdb, chain, hotspots, length, seed, helicity,
               models, advanced, paths, failure_csv):
        assert length == 65
        shutil.copy2(folder / 'trajectory-reference.pdb', Path(paths['Trajectory']) / (design_name + '.pdb'))
        print('WORKFLOW CHECK: replaying the published PDL1 trajectory boundary.', flush=True)
        return SimpleNamespace(_tmp={'best': {'aux': {'log': metrics}}},
            aux={'log': {'terminate': ''}}, get_seq=lambda **kwargs: [row['Sequence']])
    functions.binder_hallucination = replay
    sys.argv = [str(ROOT / 'bindcraft.py'), *forwarded]
    runpy.run_path(str(ROOT / 'bindcraft.py'), run_name='__main__')

def main():
    if len(sys.argv) > 1 and sys.argv[1] == '--replay-main':
        replay_main(Path(sys.argv[2]), sys.argv[3:])
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=['normal', 'accepted', 'ranking'], required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    folder = args.output_dir.resolve() / args.case
    folder.mkdir(parents=True, exist_ok=False)
    report = {'normal': normal, 'accepted': accepted, 'ranking': ranking}[args.case](folder)
    (folder / 'report.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))

if __name__ == '__main__':
    main()
