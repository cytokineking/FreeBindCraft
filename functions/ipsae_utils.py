"""
ipSAE scoring utilities for BindCraft

Based on: https://www.biorxiv.org/content/10.1101/2025.02.10.637595v1
Original implementation by Roland Dunbrack, Fox Chase Cancer Center

This module provides functions to calculate ipSAE (interface predicted 
Structural Alignment Error) scores from ColabDesign PAE matrices.

MIT license: script can be modified and redistributed for non-commercial 
and commercial use, as long as this information is reproduced.
"""

import numpy as np


def ptm_func(x, d0):
    """
    Calculate the PTM-style score for PAE values.
    
    This is the core transformation used in ipSAE/ipTM calculations.
    
    Args:
        x: PAE value(s) - can be scalar or numpy array
        d0: Normalization constant based on alignment length
    
    Returns:
        PTM-style score in range [0, 1]
    """
    return 1.0 / (1.0 + (x / d0) ** 2.0)


def calc_d0(L, pair_type='protein'):
    """
    Calculate d0 based on alignment length.
    
    From Yang and Skolnick, PROTEINS: Structure, Function, and 
    Bioinformatics 57:702–710 (2004)
    
    Args:
        L: Number of residues in the alignment
        pair_type: 'protein' or 'nucleic_acid'
    
    Returns:
        d0 value (minimum 1.0 for protein, 2.0 for nucleic acid)
    """
    L = float(max(L, 27))
    min_value = 2.0 if pair_type == 'nucleic_acid' else 1.0
    d0 = 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8
    return max(min_value, d0)


def calc_d0_array(L_array, pair_type='protein'):
    """
    Vectorized version of calc_d0 for numpy arrays.
    
    Args:
        L_array: Array of alignment lengths
        pair_type: 'protein' or 'nucleic_acid'
    
    Returns:
        Array of d0 values
    """
    L = np.array(L_array, dtype=float)
    L = np.maximum(27, L)
    min_value = 2.0 if pair_type == 'nucleic_acid' else 1.0
    return np.maximum(min_value, 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8)


def calculate_ipsae(pae_matrix, target_len, binder_len, pae_cutoff=10.0):
    """
    Calculate ipSAE score from ColabDesign PAE matrix.
    
    This implements the ipSAE_d0res calculation from the Dunbrack paper,
    which uses adaptive d0 based on the number of residues with good PAE
    values for each aligned residue.
    
    Args:
        pae_matrix: numpy array of shape (L, L) containing PAE values
                   where L = target_len + binder_len
                   Target residues are indices [0:target_len]
                   Binder residues are indices [target_len:target_len+binder_len]
        target_len: number of residues in target chain
        binder_len: number of residues in binder chain  
        pae_cutoff: PAE cutoff for considering residue pairs (default 10.0 Å)
    
    Returns:
        dict with ipSAE metrics:
            - 'ipSAE': max of binder and target direction scores (primary metric)
            - 'ipSAE_binder': max per-residue score from binder → target
            - 'ipSAE_target': max per-residue score from target → binder
            - 'n0dom': number of residues with good PAE values
    """
    total_len = target_len + binder_len
    
    # Validate input
    if pae_matrix.shape[0] != total_len or pae_matrix.shape[1] != total_len:
        raise ValueError(f"PAE matrix shape {pae_matrix.shape} does not match "
                        f"target_len ({target_len}) + binder_len ({binder_len}) = {total_len}")
    
    # Interface PAE: binder rows → target columns
    interface_pae = pae_matrix[target_len:, :target_len]
    
    # Apply PAE cutoff mask
    valid_mask = interface_pae < pae_cutoff
    
    # Count residues with good PAE values (for n0dom)
    binder_good_residues = np.any(valid_mask, axis=1).sum()
    target_good_residues = np.any(valid_mask, axis=0).sum()
    n0dom = int(binder_good_residues + target_good_residues)
    
    # Calculate per-binder-residue ipSAE scores (ipSAE_d0res)
    # For each binder residue, d0 is based on number of target residues with good PAE
    ipsae_byres = []
    for i in range(binder_len):
        valid = valid_mask[i]
        if valid.any():
            n0res = valid.sum()
            d0res = calc_d0(n0res)
            ptm_vals = ptm_func(interface_pae[i][valid], d0res)
            ipsae_byres.append(ptm_vals.mean())
        else:
            ipsae_byres.append(0.0)
    
    ipsae_byres = np.array(ipsae_byres)
    ipsae_binder_max = float(ipsae_byres.max()) if len(ipsae_byres) > 0 else 0.0
    
    # Calculate reverse direction: target rows → binder columns
    interface_pae_rev = pae_matrix[:target_len, target_len:]
    valid_mask_rev = interface_pae_rev < pae_cutoff
    
    ipsae_byres_rev = []
    for i in range(target_len):
        valid = valid_mask_rev[i]
        if valid.any():
            n0res = valid.sum()
            d0res = calc_d0(n0res)
            ptm_vals = ptm_func(interface_pae_rev[i][valid], d0res)
            ipsae_byres_rev.append(ptm_vals.mean())
        else:
            ipsae_byres_rev.append(0.0)
    
    ipsae_byres_rev = np.array(ipsae_byres_rev)
    ipsae_target_max = float(ipsae_byres_rev.max()) if len(ipsae_byres_rev) > 0 else 0.0
    
    # Take max of both directions (as in original ipSAE paper)
    ipsae = max(ipsae_binder_max, ipsae_target_max)
    
    return {
        'ipSAE': round(ipsae, 4),
        'ipSAE_binder': round(ipsae_binder_max, 4),
        'ipSAE_target': round(ipsae_target_max, 4),
        'n0dom': n0dom,
    }

