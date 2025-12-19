import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import tempfile
import os
import numpy as np
import plotly.express as px



############# HELPER FUNCTIONS #############

def extract_contact_residues(pairs):
    """Extract unique residue numbers from a list of contact pairs [( (A,res1),(B,res2) ), ... ]"""
    if not pairs or pairs == "nan":
        return []
    residues = set()
    try:
        for (a_chain, a_res), (b_chain, b_res) in pairs:
            residues.add((a_chain, int(a_res)))
            residues.add((b_chain, int(b_res)))
    except Exception:
        return []
    return list(residues)

def parse_pairs(raw_pairs, target_chain='A', binder_chain='B'):
    """
    Parse BindCraft style contact pairs.
    Example input:
    [((75,'ARG'), (28,'PHE')), ((78,'THR'), (35,'TRP'))]
    Returns:
       [('A',75), ('B',28), ('A',78), ('B',35)]
    """
    import ast

    if raw_pairs is None:
        return []

    # Convert string → Python list
    if isinstance(raw_pairs, str):
        try:
            pairs = ast.literal_eval(raw_pairs)
        except Exception:
            return []
    else:
        pairs = raw_pairs

    residues = set()

    #[((75,'ARG'), (28,'PHE')), ((78,'THR'), (35,'TRP'))]
    for left, right in pairs:
        try:
            # left = (75, 'ARG') → residue number = left[0]
            res_target = int(left[0])
            residues.add((target_chain, res_target))

            # right = (28, 'PHE') → residue number = right[0]
            res_binder = int(right[0])
            residues.add((binder_chain, res_binder))

        except Exception:
            continue

    return list(residues)

import re

def extract_design_id(pdb_filename):

    """
    Works for:
    - b_SLOG_l50_s95426_mpnn9_model2.pdb
    - 5_b_SLOG_l50_s95426_mpnn9_model2.pdb
    - 12_b_SLOG_l50_s95426_mpnn9_model1.pdb
    """
    name = os.path.splitext(os.path.basename(pdb_filename))[0]

    # remove _modelX suffix
    #name = re.sub(r'_model\d+$', '', name)
    name = "_".join(name.split("_")[:-1])

    # remove leading rank (digits + underscore)
    name = re.sub(r'^\d+_', '', name)

    return name


#######################

def highlight_residues(view, residue_list, sphere=False, cartoon_color='red', sphere_radius=1.2):
    """Highlight residues on the 3D structure as cartoon only."""
    for chain, resi in residue_list:
        try:
            view.addStyle({'chain': chain, 'resi': str(resi)}, {'cartoon': {'color': cartoon_color}})
        except Exception:
            continue
