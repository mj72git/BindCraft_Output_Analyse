import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import tempfile
import os
import numpy as np
import plotly.express as px

##### MODIFIED #####
if "analysis_done" not in st.session_state:
    st.session_state.analysis_done = False
####################
try:
    import py3Dmol
    py3dmol_available = True
except ImportError:
    py3dmol_available = False

from bindcraft_analysis_pipeline import analyze_design

st.set_page_config(page_title="BindCraft Analysis", layout="wide")
st.title("BindCraft Output Analysis (BOA) Web Application")

def format_pairs(pairs):
    if not pairs or pairs == 'nan':
        return "None"
    try:
        return " ".join([f"{a[0]} {a[1]} ↔ {b[0]} {b[1]}" for a, b in pairs])
    except Exception:
        return str(pairs)

# session_state defaults
if 'df_out' not in st.session_state:
    st.session_state.df_out = None
if 'df_rank' not in st.session_state:
    st.session_state.df_rank = None
if 'pdb_map' not in st.session_state:
    st.session_state.pdb_map = {}

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



#######################

def highlight_residues(view, residue_list, sphere=False, cartoon_color='red', sphere_radius=1.2):
    """Highlight residues on the 3D structure as cartoon only."""
    for chain, resi in residue_list:
        try:
            view.addStyle({'chain': chain, 'resi': str(resi)}, {'cartoon': {'color': cartoon_color}})
        except Exception:
            continue

############# SIDEBAR SETTINGS #############
st.sidebar.header("Settings")
target_chain = st.sidebar.text_input("Target Chain Letter", value="A")
binder_chain = st.sidebar.text_input("Binder Chain Letter", value="B")

add_target_res_offset = st.sidebar.number_input("Target Residue Offset ", value=0)
st.sidebar.markdown("##### (e.g. your Output BindCraft target chain starts from Residue 1 but your initial target chain, starts from 24. you should enter 23)")
blankk = st.sidebar.header("")

mj = st.sidebar.header("App created by MJ Shadfar")
st.sidebar.caption("BOA v1.2")

############# FILE UPLOAD & ANALYSIS #############
if not st.session_state.analysis_done:
    st.subheader("Upload Files")
    pdb_files = st.file_uploader("Upload PDB files", type=["pdb", "cif"], accept_multiple_files=True)
    csv_file = st.file_uploader("Upload design metrics CSV. (final_design_stats.csv)", type="csv")

    if pdb_files and csv_file:
        st.success("Files uploaded successfully.")
        if st.button("Run Analysis"):
            tmpdir = tempfile.mkdtemp()
            freesasa_available = (os.system("which freesasa > /dev/null") == 0)
            csv_path = os.path.join(tmpdir, "metrics.csv")
            with open(csv_path, "wb") as f:
                f.write(csv_file.read())
            df_metrics = pd.read_csv(csv_path)

            results = []
            progress = st.progress(0)
            for i, pdb in enumerate(pdb_files):
                pdb_path = os.path.join(tmpdir, pdb.name)
                with open(pdb_path, "wb") as f:
                    f.write(pdb.read())

                r = analyze_design(pdb_path, target_chain=target_chain, binder_chain=binder_chain,
                                   add_target_res_offset=add_target_res_offset,
                                   freesasa_available=freesasa_available, tmpdir=tmpdir)

                base = os.path.splitext(os.path.basename(pdb_path))[0]
                base = "_".join(base.split("_")[:-1]) or os.path.splitext(os.path.basename(pdb_path))[0]

                matched_row = {}
                if 'Design' in df_metrics.columns:
                    m = df_metrics[df_metrics['Design'].astype(str).str.contains(base)]
                    if len(m) == 1:
                        matched_row = m.iloc[0].to_dict()
                    else:
                        m = df_metrics[df_metrics['Design'].astype(str) == base]
                        if len(m) == 1:
                            matched_row = m.iloc[0].to_dict()

                try:
                    with open(pdb_path, "r") as fh:
                        pdb_text = fh.read()
                except UnicodeDecodeError:
                    with open(pdb_path, "r", encoding="latin-1") as fh:
                        pdb_text = fh.read()
                st.session_state.pdb_map[base] = pdb_text

                record = {'design_id': base}
                for col in ['Average_pLDDT','Average_i_pLDDT','Average_pTM','Average_i_pAE','Average_i_pTM','Average_pAE','Average_dG','Average_dSASA','Average_Binder_pLDDT','Average_n_InterfaceResidues']:
                    record[col] = matched_row.get(col, np.nan)
                record.update({
                    'n_contacts_3A': r.get('n_contacts_3A'),
                    'n_contacts_4A': r.get('n_contacts_4A'),
                    'n_target_interface_residues': r.get('n_target_interface_residues'),
                    'n_binder_interface_residues': r.get('n_binder_interface_residues'),
                    'hbond_like_count': r.get('hbond_like_count'),
                    'clash_count': r.get('clash_count'),
                    'dsasa': r.get('dsasa'),
                    'target_seq': r.get('target_seq'),
                    'binder_seq': r.get('binder_seq'),
                    'pairs_3A': r.get('pairs_3A'),
                    'pairs_4A': r.get('pairs_4A'),
                    'hbond_pairs': r.get('hbond_pairs')
                })
                results.append(record)
                progress.progress((i+1)/len(pdb_files))

            st.session_state.df_out = pd.DataFrame(results)
            df_rank = st.session_state.df_out.copy()
            for col in ['Average_i_pTM','dsasa','Average_pLDDT']:
                if col not in df_rank.columns:
                    df_rank[col] = np.nan
            st.session_state.df_rank = df_rank.sort_values(by=['Average_i_pTM','dsasa','Average_pLDDT'], ascending=[False, False, False])
            st.session_state.analysis_done = True
            st.rerun()

############# AFTER ANALYSIS #############
if st.session_state.analysis_done and st.session_state.df_out is not None:
    tab1, tab2, tab3 = st.tabs(["Summary", "Visualizations", "Details"])

    with tab1:
        st.subheader("Analysis Summary")
        df_display = st.session_state.df_out.copy()
        for col in ['pairs_3A', 'pairs_4A', 'hbond_pairs']:
            if col in df_display.columns:
                df_display[col] = df_display[col].apply(lambda x: str(x))
        st.dataframe(df_display)

        st.subheader("Ranked Designs")
        df_rank_display = st.session_state.df_rank.copy()
        for col in ['pairs_3A', 'pairs_4A', 'hbond_pairs']:
            if col in df_rank_display.columns:
                df_rank_display[col] = df_rank_display[col].apply(lambda x: str(x))
        st.dataframe(df_rank_display)

        st.download_button("Download Summary CSV", st.session_state.df_out.to_csv(index=False), "bindcraft_analysis_summary.csv")
        st.download_button("Download Ranked CSV", st.session_state.df_rank.to_csv(index=False), "bindcraft_ranked.csv")

        if st.button("Reset"):
            st.session_state.df_out = None
            st.session_state.df_rank = None
            st.session_state.pdb_map = {}
            st.session_state.analysis_done = False
            st.rerun()

    with tab2:
        st.subheader("Visualizations")

        # ensure df is defined
        df = st.session_state.df_out

        # -------------------------------
        # BAR PLOT (4Å contacts)
        # -------------------------------
        if df is not None and df.shape[0] > 0:
            fig_contacts = px.bar(
                df,
                x='design_id',
                y='n_contacts_4A',
                title='Number of Contacts (4Å)'
            )
            st.plotly_chart(fig_contacts, use_container_width=True)
        else:
            st.info("No data to plot. Please run analysis first.")

        # -------------------------------
        # 3D STRUCTURE VIEWER
        # -------------------------------
        st.subheader("3D Structure Viewer")

        if not py3dmol_available:
            st.error("py3Dmol is not installed.")
        else:
            if df is None or df.shape[0] == 0:
                st.info("No designs available. Please run analysis first.")
            else:
                # Design selection
                selected_design = st.selectbox("Select a design", df['design_id'].tolist())
                row = df[df['design_id'] == selected_design].iloc[0]

                # Handle raw_pairs
                raw_pairs = row.get('pairs_4A')
                if isinstance(raw_pairs, str) and raw_pairs.strip().lower() in ["nan", "", "none"]:
                    raw_pairs = None

                # Parse contacts
                contacts = parse_pairs(raw_pairs, target_chain=target_chain, binder_chain=binder_chain)

                # Apply target offset and normalize chains
                contacts_offset = []
                t_chain = str(target_chain).strip().upper()
                b_chain = str(binder_chain).strip().upper()

                for c, r in contacts:
                    try:
                        chain_norm = str(c).strip().upper()
                        resi = int(r)
                        if chain_norm == t_chain:
                            resi += int(add_target_res_offset)
                        contacts_offset.append((chain_norm, resi))
                    except:
                        continue

                # Prepare highlight lists for target and binder
                highlight_target = [(c, r) for c, r in contacts_offset if c == t_chain]
                highlight_binder = [(c, r) for c, r in contacts_offset if c == b_chain]

                pdb_text = st.session_state.pdb_map[selected_design]

                # Chain visibility + colors
                show_target = st.checkbox("Show Target Chain", value=True, key=f"show_target_{selected_design}")
                show_binder = st.checkbox("Show Binder Chain", value=True, key=f"show_binder_{selected_design}")
                highlight_binding = st.checkbox("Highlight binding spot", value=True,
                                                key=f"highlight_{selected_design}")

                target_color = st.selectbox("Target Color", ["limegreen", "red", "orange", "magenta", "yellow", "cyan"],
                                            index=0, key=f"tcolor_{selected_design}")
                binder_color = st.selectbox("Binder Color",
                                            ["deepskyblue", "red", "orange", "magenta", "yellow", "cyan"], index=0,
                                            key=f"bcolor_{selected_design}")

                # Create viewer
                view = py3Dmol.view(width=900, height=500)
                view.addModel(pdb_text, 'pdb')
                view.setBackgroundColor('0x303030')

                # Cartoon visibility control
                if show_target:
                    view.setStyle({'chain': t_chain}, {'cartoon': {'color': target_color}})
                else:
                    view.setStyle({'chain': t_chain}, {'cartoon': {'opacity': 0.0}})

                if show_binder:
                    view.setStyle({'chain': b_chain}, {'cartoon': {'color': binder_color}})
                else:
                    view.setStyle({'chain': b_chain}, {'cartoon': {'opacity': 0.0}})

                # Contact highlighting, only for visible chains and if checkbox is checked
                if highlight_binding:
                    if show_target:
                        highlight_residues(
                            view,
                            highlight_target,
                            sphere=True,
                            cartoon_color='purple',
                            sphere_radius=1.3
                        )

                    if show_binder:
                        highlight_residues(
                            view,
                            highlight_binder,
                            sphere=True,
                            cartoon_color='purple',
                            sphere_radius=1.3
                        )

                view.zoomTo()
                components.html(view._make_html(), height=500, width=900)

    with tab3:
        st.subheader("Per-Design Details")
        for _, row in st.session_state.df_out.iterrows():
            with st.expander(f"Details for {row['design_id']}"):
                st.write("**Sequences**")
                st.text(f"Target: {row['target_seq']}")
                st.text(f"Binder: {row['binder_seq']}")
                st.write("**Contacts (3Å)**")
                st.text(format_pairs(row['pairs_3A']))
                st.write("**Contacts (4Å)**")
                st.text(format_pairs(row['pairs_4A']))
                st.write("**H-bond-like pairs**")
                st.text(format_pairs(row['hbond_pairs']))
else:
    st.info("Please upload PDB files and a CSV to start analysis.")
    blankk = st.header("")
    st.markdown("#### Before starting, it is better to download and read the tutorial")
    file_path = "help.txt"
    with open(file_path, "r") as ff:
        txt_content = ff.read()
    st.download_button(
        label="Download HELP File",
        data=txt_content,
        file_name="help.txt",
        mime="text/plain"
    )
