#!/usr/bin/env python3
"""Exercise native scoring/relaxation and GPU MPNN on a bundled PDL1 benchmark.

This uses an existing published design, not a newly generated candidate. Run:
JAX_PLATFORMS=metal .venv-metaljax/bin/python extras/check_apple_installation.py
"""
import os
import sys
from pathlib import Path
root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root))
os.chdir(root)
os.environ.setdefault('JAX_PLATFORMS', 'metal')
os.environ['OPENMM_PLATFORM_ORDER'] = 'OpenCL'
os.environ['LIBCIFPP_DATA_DIR'] = str(root/'.local/dssp/share/libcifpp')
os.environ['FASPR_BIN'] = str(root/'.local/src/FASPR/FASPR')
os.environ['SC_RS_BIN'] = str(root/'.local/src/sc-rs/target/release/sc')
import json
import logging
import zipfile
from Bio.PDB import PDBParser, DSSP
from functions import check_jax_gpu, calc_ss_percentage, score_interface, mpnn_gen_sequence
from functions.pr_alternative_utils import openmm_relax

logging.basicConfig(level=logging.WARNING)
logging.getLogger('functions').setLevel(logging.DEBUG)
check_jax_gpu()
folder = root/'runs/native-validation'
folder.mkdir(parents=True, exist_ok=True)
pdb = str(folder/'reference.pdb')
out = str(folder/'relaxed.pdb')
with zipfile.ZipFile('performance_data/pdl1_miniprotein/PDL1-freebindcraft.zip') as archive:
    Path(pdb).write_bytes(archive.read('8_PDL1_l65_s321660_mpnn17_model1.pdb'))
a = json.loads((root/'settings_advanced/default_4stage_multimer.json').read_text())
a.update(dssp_path=str(root/'.local/dssp/bin/mkdssp'), num_seqs=1, max_mpnn_sequences=1)
dssp = DSSP(PDBParser(QUIET=True).get_structure('p',pdb)[0], pdb, dssp=a['dssp_path'])
assert len(dssp) == 180, len(dssp)
ss = calc_ss_percentage(pdb, a, 'B')
assert sum(ss[:3]) > 99
perf = {}
plat = openmm_relax(pdb, out, override_platform_order=['OpenCL'], perf_report=perf)
assert plat == 'OpenCL' and perf['faspr_success'] and perf['post_faspr_min_platform'] == 'OpenCL', perf
metrics, aa, residues = score_interface(out, 'B', use_pyrosetta=False)
assert metrics['interface_nres'] > 0 and residues, metrics
seqs = mpnn_gen_sequence(out, 'B', residues, a)
assert len(seqs['seq']) == 1 and len(seqs['seq'][0].replace('/', '')) == 180
report = dict(reference='Bundled benchmark complex, not the new design',
              dssp_residues=len(dssp), secondary_structure=ss[:3], relaxation=perf,
              interface=metrics, mpnn_sequence=seqs['seq'][0][-65:])
(folder/'report.json').write_text(json.dumps(report, indent=2, default=float)+'\n')
print(json.dumps(report, indent=2, default=float))
