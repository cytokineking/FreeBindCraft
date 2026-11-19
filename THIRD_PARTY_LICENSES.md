# Third-Party Licenses and Upstream Notices

FreeBindCraft incorporates or builds upon upstream software, bundled tools, and runtime libraries under their respective licenses, preserved below.

---

## 1. BindCraft (Upstream Repository)

Source: https://github.com/martinpacesa/BindCraft  
Original notice from upstream `LICENSE`:

```text
MIT License

Copyright (c) 2024 Martin Pacesa

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

## 2. ipSAE (`ipSAE.py`, `functions/ipsae_utils.py`)

Adapted from Roland Dunbrack (Fox Chase Cancer Center):  
https://www.biorxiv.org/content/10.1101/2025.02.10.637595v1

Original header notice (`ipSAE.py` / `functions/ipsae_utils.py`):

```text
Roland Dunbrack
Fox Chase Cancer Center
version 3
April 6, 2025
MIT license: script can be modified and redistributed for non-commercial and commercial use, as long as this information is reproduced.
```

---

## 3. FASPR (`FASPR-master/`, `functions/FASPR`)

Source: https://github.com/tommyhuangthu/FASPR  
Verbatim notice from `FASPR-master/LICENSE`:

```text
MIT License

Copyright (c) 2020 Xiaoqiang Huang

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

## 4. ColabDesign (`ColabDesign-main/`)

Source: https://github.com/sokrypton/ColabDesign  
Verbatim notice from `ColabDesign-main/LICENSE.txt`:

```text
/*
 * ----------------------------------------------------------------------------
 * "THE BEER-WARE LICENSE" (Revision 42):
 * <so3@mit.edu> wrote this file.  As long as you retain this notice you
 * can do whatever you want with this stuff. If we meet some day, and you think
 * this stuff is worth it, you can buy me a beer in return.  Sergey Ovchinnikov
 * ----------------------------------------------------------------------------
 */
```

---

## 5. Additional External Tools and Runtime Dependencies

| Component | Repository / Source | License |
| --- | --- | --- |
| `tmol` & `tmol-ops-cu12` | https://github.com/RosettaCommons/tmol | BSD-3-Clause (`Copyright (c) 2026, Institute for Protein Design, University of Washington`) |
| `sc-rs` (`functions/sc`) | https://github.com/cytokineking/sc-rs | MIT |
| `FreeSASA` | https://github.com/mittinatten/freesasa | MIT |
| `PDBFixer` | https://github.com/openmm/pdbfixer | MIT |
| `PyTorch` | https://pytorch.org | BSD-3-Clause |
| `JAX` | https://github.com/jax-ml/jax | Apache-2.0 |
| `OpenMM` | https://openmm.org | Open-source (PSF-like / MIT / LGPL per upstream distribution) |
| `Biopython` | https://biopython.org | Biopython License Agreement |
| `DSSP` (`functions/dssp`) | https://github.com/PDB-REDO/dssp | Upstream license for bundled standalone binary |
| `PyRosetta` & `functions/DAlphaBall.gcc` *(optional)* | https://www.pyrosetta.org | Rosetta Commons License (separate commercial license required for non-academic use) |
