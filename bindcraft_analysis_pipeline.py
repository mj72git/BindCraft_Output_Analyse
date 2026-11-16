import os
import sys
import argparse
import subprocess
from shutil import which
import math
import textwrap

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib import distances



# ----------------------------- Analysis helpers -----------------------------


def residue_contacts(target_atoms, binder_atoms, cutoff=4.0, add_target_res_offset=0):
    """
    Returns:
    - contacts: sorted list of ((resA_id, resA_name), (resB_id, resB_name))
    - hydrophobic_patches: list of contacts specifically between two hydrophobic residues
    """
    # 1. Define Hydrophobic Lookup Table
    HYDROPHOBIC_AAS = {'VAL', 'ILE', 'LEU', 'PHE', 'MET', 'TRP', 'ALA', 'PRO'}

    if len(target_atoms) == 0 or len(binder_atoms) == 0:
        return [], []

    # 2. Calculate Distances
    d = distances.distance_array(target_atoms.positions, binder_atoms.positions)
    i, j = np.where(d <= cutoff)
    
    pairs = set()
    hydrophobic_patches = []

    # 3. Analyze Contacts
    for a, b in zip(i, j):
        # Extract residue info
        resA_id = int(target_atoms[a].resid + add_target_res_offset)
        resA_name = target_atoms[a].resname
        
        resB_id = int(binder_atoms[b].resid)
        resB_name = binder_atoms[b].resname
        
        resA_info = (resA_id, resA_name)
        resB_info = (resB_id, resB_name)
        
        # Add to unique set of pairs
        if (resA_info, resB_info) not in pairs:
            pairs.add((resA_info, resB_info))
            
            # 4. Identify Stabilizing Patches
            if resA_name in HYDROPHOBIC_AAS and resB_name in HYDROPHOBIC_AAS:
                hydrophobic_patches.append((resA_info, resB_info))

    sorted_contacts = sorted(list(pairs), key=lambda x: (x[0][0], x[1][0]))
    sorted_patches = sorted(hydrophobic_patches, key=lambda x: (x[0][0], x[1][0]))

    return sorted_contacts, sorted_patches


def count_hbond_like(target_atoms, binder_atoms, cutoff=3.5):
    """Count N/O atom pairs within cutoff. Return (count, unique_pairs_list).
       pairs reported as (target_resid, binder_resid)
    """
    targ_NO = target_atoms.select_atoms("name N O")
    bind_NO = binder_atoms.select_atoms("name N O")
    if len(targ_NO) == 0 or len(bind_NO) == 0:
        return 0, []
    d = distances.distance_array(targ_NO.positions, bind_NO.positions)
    i, j = np.where(d <= cutoff)
    pairs = set()
    for a, b in zip(i, j):
        pairs.add((int(targ_NO[a].resid), int(bind_NO[b].resid)))
    return len(pairs), sorted(pairs)





def extract_sequence(atomgroup):
    """Return one-letter sequence for an MDAnalysis atomgroup's residues (unique residues in order)."""
    # MDAnalysis resid/resname mapping
    three_to_one = {
        'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
        'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
        'THR':'T','TRP':'W','TYR':'Y','VAL':'V'
    }
    seq = []
    seen = set()
    for res in atomgroup.residues:
        if res.resid in seen:
            continue
        seen.add(res.resid)
        aa = three_to_one.get(res.resname.strip().upper(), 'X')
        seq.append(aa)
    return ''.join(seq)


# ----------------------------- Main design analysis -----------------------------

def analyze_design(pdb_path, target_chain='A', binder_chain='B', add_target_res_offset=0, tmpdir=None):
    """Analyze a single design PDB file. Returns dict of computed metrics."""
    out = {}
    out['pdb'] = pdb_path
    try:
        u = mda.Universe(pdb_path)
    except Exception as e:
        out['error'] = f"MDAnalysis load error: {e}"
        return out

    # Support chain selection: MDAnalysis may use chainIDs as segids depending on file
    # We'll try both selectors (chainID and segid) for robustness
    def select_non_h(target_chain_letter):
        sel1 = u.select_atoms(f"chainID {target_chain_letter} and not name H*")
        if len(sel1) > 0:
            return sel1
        sel2 = u.select_atoms(f"segid {target_chain_letter} and not name H*")
        if len(sel2) > 0:
            return sel2
        # fallback: try chain by residue range? return empty
        return sel1

    target = select_non_h(target_chain)
    binder = select_non_h(binder_chain)

    out['target_atoms'] = len(target)
    out['binder_atoms'] = len(binder)

    # Contacts
   # 3.0A Contacts and Hydrophobic Patches
    pairs3, hypho_list = residue_contacts(target, binder, cutoff=3.0, add_target_res_offset=add_target_res_offset)
    
    # 4.0A Contacts
    pairs4, _ = residue_contacts(target, binder, cutoff=4.0, add_target_res_offset=add_target_res_offset)
    out['n_contacts_3A'] = len(pairs3)
    out['n_contacts_4A'] = len(pairs4)
    out['pairs_3A'] = pairs3
    out['pairs_4A'] = pairs4
    out['hypho'] = hypho_list

    # Unique contacting residues
    target_res_set = set([p[0][0] for p in pairs3])
    binder_res_set = set([p[1][0] for p in pairs3])
    out['n_target_interface_residues'] = len(target_res_set)
    out['n_binder_interface_residues'] = len(binder_res_set)


    # Sequences
    out['target_seq'] = extract_sequence(target)
    out['binder_seq'] = extract_sequence(binder)

    return out

