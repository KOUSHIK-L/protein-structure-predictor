import streamlit as st
import py3Dmol
import requests
import biotite.structure as struc
from biotite.structure.io.pdb import PDBFile
from streamlit_option_menu import option_menu
import re
import numpy as np
import pandas as pd
from io import StringIO

# Page config
st.set_page_config(
    page_title="ESMFold | Protein Structure Predictor",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# CSS
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');

html, body, [class*="css"] { font-family: 'Inter', sans-serif !important; }

/* Hide default Streamlit chrome */
#MainMenu { visibility: hidden; }
footer    { visibility: hidden; }
header [data-testid="stToolbar"] { visibility: hidden; }

/* Sidebar */
[data-testid="stSidebar"] {
    background: #ffffff;
    border-right: 1px solid #e2e8f0;
}

/* Predict button */
[data-testid="stSidebar"] .stButton > button {
    background: linear-gradient(135deg, #0f3460 0%, #1a6eb5 100%) !important;
    color: white !important;
    border: none !important;
    border-radius: 10px !important;
    padding: 0.7rem 1rem !important;
    font-weight: 600 !important;
    font-size: 0.95rem !important;
    width: 100% !important;
    box-shadow: 0 2px 10px rgba(15,52,96,0.3) !important;
    transition: all 0.2s ease !important;
}
[data-testid="stSidebar"] .stButton > button:hover {
    opacity: 0.88 !important;
    box-shadow: 0 4px 18px rgba(15,52,96,0.4) !important;
    transform: translateY(-1px) !important;
}
[data-testid="stSidebar"] .stButton > button:disabled {
    background: #e2e8f0 !important;
    color: #94a3b8 !important;
    box-shadow: none !important;
    transform: none !important;
}

/* Download button */
.stDownloadButton > button {
    background: #f0f9ff !important;
    color: #0f3460 !important;
    border: 1px solid #bae6fd !important;
    border-radius: 8px !important;
    font-weight: 500 !important;
    transition: background 0.2s !important;
}
.stDownloadButton > button:hover { background: #e0f2fe !important; }

/* Tabs */
button[role="tab"] { font-weight: 500 !important; }
button[role="tab"][aria-selected="true"] {
    color: #0f3460 !important;
    font-weight: 600 !important;
}
</style>
""", unsafe_allow_html=True)

# Session state
_defaults = {"txt": None, "uniprot_meta": None, "pdb_string": None, "prediction_error": None}
for k, v in _defaults.items():
    if k not in st.session_state:
        st.session_state[k] = v


# HELPERS

def render_mol(pdb, color_mode="Spectrum"):
    pdbview = py3Dmol.view(width=800, height=520)
    pdbview.addModel(pdb, "pdb")
    if color_mode == "Spectrum":
        pdbview.setStyle({"cartoon": {"color": "spectrum"}})
    elif color_mode == "plDDT Confidence":
        pdbview.setStyle({"cartoon": {
            "colorscheme": {"prop": "b", "gradient": "rwb", "min": 50, "max": 100}
        }})
    elif color_mode == "Secondary Structure":
        pdbview.setStyle({"cartoon": {"color": "ssJmol"}})
    pdbview.setBackgroundColor("white")
    pdbview.zoomTo()
    pdbview.zoom(2, 800)
    pdbview.spin(True)
    st.iframe(pdbview._make_html(), height=520)


def is_valid_protein_sequence(seq):
    return bool(re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", seq))


def parse_fasta(text):
    return "".join(
        line.strip() for line in text.splitlines()
        if not line.strip().startswith(">")
    )


def plddt_badge(score):
    if score >= 90:
        bg, fg, label = "#dbeafe", "#1d4ed8", "Very High"
    elif score >= 70:
        bg, fg, label = "#d1fae5", "#065f46", "Confident"
    elif score >= 50:
        bg, fg, label = "#fef9c3", "#854d0e", "Low"
    else:
        bg, fg, label = "#fee2e2", "#991b1b", "Very Low"
    return (
        f'<span style="background:{bg};color:{fg};padding:0.3rem 0.9rem;'
        f'border-radius:999px;font-weight:600;font-size:0.88rem;'
        f'display:inline-block">{label} &nbsp;·&nbsp; {score:.1f}</span>'
    )


def metric_card(label, value, sub=""):
    return f"""
<div style="background:#f8fafc;border:1px solid #e2e8f0;border-radius:12px;
            padding:1.2rem;text-align:center;height:100%">
    <div style="font-size:0.7rem;font-weight:600;color:#64748b;
                text-transform:uppercase;letter-spacing:0.8px">{label}</div>
    <div style="font-size:1.9rem;font-weight:700;color:#0f3460;
                margin:0.3rem 0 0.1rem;line-height:1.1">{value}</div>
    <div style="font-size:0.78rem;color:#94a3b8">{sub}</div>
</div>"""


def gradient_header(title, subtitle=""):
    sub_html = (
        f'<p style="margin:0.4rem 0 0;opacity:0.8;font-size:0.88rem">{subtitle}</p>'
        if subtitle else ""
    )
    return f"""
<div style="background:linear-gradient(135deg,#0f3460 0%,#1a6eb5 100%);
            border-radius:16px;padding:2rem 2.5rem;color:white;margin-bottom:1.5rem">
    <h2 style="margin:0;font-size:1.55rem;font-weight:700">{title}</h2>
    {sub_html}
</div>"""


def fetch_fasta_from_url(url):
    try:
        r = requests.get(url, timeout=10)
        return (r.text, None) if r.status_code == 200 else (None, f"HTTP {r.status_code}")
    except Exception as e:
        return None, str(e)


def extract_uniprot_accession(url):
    m = re.search(r"[/=]([A-Z][0-9][A-Z0-9]{3}[0-9])\b", url, re.IGNORECASE)
    return m.group(1).upper() if m else None


def fetch_uniprot_metadata(accession):
    try:
        r = requests.get(f"https://rest.uniprot.org/uniprotkb/{accession}.json", timeout=10)
        if r.status_code != 200:
            return None
        d = r.json()
        name = (d.get("proteinDescription", {}).get("recommendedName", {})
                 .get("fullName", {}).get("value", "Unknown"))
        organism = d.get("organism", {}).get("scientificName", "Unknown")
        gene = next((g.get("geneName", {}).get("value", "") for g in d.get("genes", [])), "")
        fn = next(
            (t.get("value", "") for c in d.get("comments", [])
             if c.get("commentType") == "FUNCTION" for t in c.get("texts", [])),
            "Not available",
        )
        return {"name": name, "gene": gene, "organism": organism, "function": fn, "accession": accession}
    except Exception:
        return None


def extract_plddt_per_residue(pdb_string):
    residues = {}
    for line in pdb_string.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                residues[int(line[22:26])] = float(line[60:66])
            except ValueError:
                pass
    if not residues:
        return pd.DataFrame()
    items = sorted(residues.items())
    return pd.DataFrame({"plDDT": [v for _, v in items]}, index=[k for k, _ in items])


def get_secondary_structure(pdb_string):
    try:
        pdb_file = PDBFile.read(StringIO(pdb_string))
        atom_array = pdb_file.get_structure(model=1)
        protein = atom_array[struc.filter_amino_acids(atom_array)]
        sse = struc.annotate_sse(protein)
        return {
            "Alpha Helix": int(np.sum(sse == "a")),
            "Beta Sheet":  int(np.sum(sse == "b")),
            "Coil / Loop": int(np.sum(sse == "c")),
        }
    except Exception:
        return None


def get_aa_composition(sequence):
    aa_names = {
        "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
        "E": "Glu", "Q": "Gln", "G": "Gly", "H": "His", "I": "Ile",
        "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
        "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    }
    counts = {aa: sequence.count(aa) for aa in aa_names}
    return pd.DataFrame(
        {"Count": counts.values()},
        index=[f"{aa} ({aa_names[aa]})" for aa in counts],
    ).sort_values("Count", ascending=False)


@st.cache_data(show_spinner=False)
def call_esmfold_api(sequence):
    try:
        r = requests.post(
            "https://api.esmatlas.com/foldSequence/v1/pdb/",
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            data=sequence,
            timeout=120,
        )
    except requests.exceptions.Timeout:
        return None, "Request timed out. Try a shorter sequence or try again later."
    except requests.exceptions.ConnectionError:
        return None, "Cannot reach the ESMFold API. It may be temporarily down."
    if r.status_code != 200:
        return None, f"ESMFold API returned HTTP {r.status_code}. Try again later."
    return r.content.decode("utf-8"), None


# SIDEBAR

with st.sidebar:
    # Brand
    st.markdown("""
    <div style="text-align:center;padding:1.2rem 0 0.8rem">
        <div style="font-size:2.8rem">🧬</div>
        <div style="font-size:1.15rem;font-weight:700;color:#0f3460;margin-top:0.2rem">
            ESMFold
        </div>
        <div style="font-size:0.75rem;color:#94a3b8;margin-top:0.1rem">
            Protein Structure Predictor
        </div>
    </div>
    """, unsafe_allow_html=True)

    st.divider()

    # ① Input
    st.markdown(
        '<p style="font-size:0.7rem;font-weight:600;color:#64748b;'
        'text-transform:uppercase;letter-spacing:0.8px;margin-bottom:0.4rem">'
        '① Input Sequence</p>',
        unsafe_allow_html=True,
    )

    txt_input = st.text_area(
        "Paste amino acid sequence",
        placeholder="ACDEFGHIKLMNPQRSTVWY...",
        height=95,
        label_visibility="collapsed",
    )
    if txt_input:
        clean = txt_input.strip().upper()
        if clean != st.session_state.txt:
            st.session_state.txt = clean
            st.session_state.uniprot_meta = None
            st.session_state.pdb_string = None

    fasta_url = st.text_input(
        "UniProt FASTA URL",
        placeholder="https://rest.uniprot.org/uniprotkb/P69905.fasta",
    )
    if fasta_url:
        fasta_text, err = fetch_fasta_from_url(fasta_url)
        if err:
            st.error(err)
        else:
            seq = parse_fasta(fasta_text).upper()
            if seq != st.session_state.txt:
                st.session_state.txt = seq
                st.session_state.pdb_string = None
                acc = extract_uniprot_accession(fasta_url)
                st.session_state.uniprot_meta = fetch_uniprot_metadata(acc) if acc else None

    uploaded = st.file_uploader("Upload FASTA / TXT", type=["fasta", "fa", "txt"])
    if uploaded:
        seq = parse_fasta(uploaded.read().decode()).upper()
        if seq != st.session_state.txt:
            st.session_state.txt = seq
            st.session_state.uniprot_meta = None
            st.session_state.pdb_string = None

    txt = st.session_state.txt

    # Validation + quick stats
    if txt:
        if not is_valid_protein_sequence(txt):
            st.error("Invalid characters - use standard single-letter AA codes.")
            st.session_state.txt = None
            txt = None
        else:
            c1, c2 = st.columns(2)
            c1.metric("Length", f"{len(txt)} aa")
            c2.metric("Unique AAs", f"{len(set(txt))}/20")
            if len(txt) > 400:
                st.warning("Sequences >400 aa may time out or reduce accuracy.")

    st.divider()

    # ② Options
    st.markdown(
        '<p style="font-size:0.7rem;font-weight:600;color:#64748b;'
        'text-transform:uppercase;letter-spacing:0.8px;margin-bottom:0.4rem">'
        '② Options</p>',
        unsafe_allow_html=True,
    )
    max_len = st.slider("Max residues to fold", 50, 600, 400,
                        help="ESMFold is most accurate at ≤400 aa")
    color_mode = st.selectbox(
        "3D color scheme",
        ["Spectrum", "plDDT Confidence", "Secondary Structure"],
        help="plDDT: blue = high confidence, red = low",
    )

    st.divider()

    # Predict
    predict = st.button("🔬  Predict Structure", disabled=(txt is None), width='stretch')

    if predict and txt:
        with st.spinner("Running ESMFold..."):
            pdb, err = call_esmfold_api(txt[:max_len])
        if err:
            st.session_state.prediction_error = err
            st.session_state.pdb_string = None
        else:
            st.session_state.pdb_string = pdb
            st.session_state.prediction_error = None

# main area section

page = option_menu(
    None,
    ["About", "Structure Prediction"],
    icons=["info-circle", "cpu"],
    orientation="horizontal",
    styles={
        "container":        {"padding": "0", "background-color": "#f8fafc"},
        "icon":             {"color": "#64748b"},
        "nav-link":         {"font-size": "0.9rem", "font-weight": "500", "color": "#64748b"},
        "nav-link-selected": {
            "background-color": "#0f3460",
            "color": "white",
            "font-weight": "600",
        },
    },
)


# ABOUT PAGE

if page == "About":
    st.markdown(gradient_header(
        "🧬 About ESMFold",
        "AI-powered 3D protein structure prediction from sequence alone"
    ), unsafe_allow_html=True)

    col1, col2 = st.columns(2)

    with col1:
        with st.container(border=True):
            st.markdown("#### ⚙️ How It Works")
            st.markdown("""
1. **ESM-2 transformer** reads the amino acid sequence and encodes evolutionary context learned from 250 million proteins
2. **Folding head** converts sequence embeddings into 3D atomic coordinates
3. **plDDT score** quantifies per-residue prediction confidence (0–100)
4. Full prediction runs in **seconds to minutes** - no MSA database needed
            """)

    with col2:
        with st.container(border=True):
            st.markdown("#### 📊 plDDT Confidence Guide")
            for bg, rng, label in [
                ("#dbeafe", "> 90",  "Very High - backbone is reliable"),
                ("#d1fae5", "70–90", "Confident - generally accurate"),
                ("#fef9c3", "50–70", "Low - treat with caution"),
                ("#fee2e2", "< 50",  "Very Low - likely disordered"),
            ]:
                st.markdown(
                    f'<div style="background:{bg};border-radius:8px;padding:0.5rem 0.85rem;'
                    f'margin-bottom:0.4rem;font-size:0.88rem"><b>{rng}</b> - {label}</div>',
                    unsafe_allow_html=True,
                )

    with st.container(border=True):
        st.markdown("#### 📥 Accepted Input Formats")
        c1, c2, c3 = st.columns(3)
        with c1:
            st.markdown("**Direct sequence**\nPaste raw single-letter amino acid codes into the sidebar")
        with c2:
            st.markdown("**UniProt FASTA URL**\nAutomatically fetches the sequence and protein metadata")
        with c3:
            st.markdown("**File upload**\nSupports `.fasta`, `.fa`, and `.txt` files")

    with st.container(border=True):
        st.markdown("#### 📚 References")
        st.markdown("""
- Lin et al. (2023). *Evolutionary-scale prediction of atomic-level protein structure with a language model.* [DOI: 10.1126/science.ade2574](https://www.science.org/doi/10.1126/science.ade2574)
- [ESM Metagenomic Atlas - Meta AI blog](https://ai.facebook.com/blog/protein-folding-esmfold-metagenomics/)
- [AlphaFold vs ESMFold - Nature news](https://www.nature.com/articles/d41586-022-03539-1)
- [Tutorial by Chanin Nantasenamat (dataprofessor) - YouTube](https://youtu.be/GHoE4VkDehY?si=aKNldfYGv9eW1I7f)
""")


# PREDICTION PAGE

if page == "Structure Prediction":

    st.markdown(gradient_header(
        "🔬 Structure Prediction",
        "Enter a sequence in the sidebar → set options → click Predict Structure",
    ), unsafe_allow_html=True)

    # Error state
    if st.session_state.prediction_error:
        st.error(st.session_state.prediction_error)
        st.info(
            "If the ESMFold API is unavailable, try "
            "[ColabFold](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb) "
            "as an alternative."
        )

    # UniProt metadata card
    if st.session_state.uniprot_meta:
        meta = st.session_state.uniprot_meta
        with st.container(border=True):
            st.markdown("**🔎 UniProt Protein Info**")
            c1, c2, c3 = st.columns(3)
            c1.markdown(f"**Protein**\n{meta['name']}")
            c2.markdown(f"**Gene**\n{meta['gene'] or '-'}")
            c3.markdown(f"**Organism**\n*{meta['organism']}*")
            if meta["function"] != "Not available":
                st.caption(f"**Function:** {meta['function']}")

    # Empty / Welcome state
    if not st.session_state.pdb_string and not st.session_state.prediction_error:
        st.markdown("""
        <div style="text-align:center;padding:2.5rem 1rem 1rem;color:#94a3b8">
            <div style="font-size:3.5rem">🔬</div>
            <div style="font-size:1.1rem;font-weight:600;color:#64748b;margin-top:0.6rem">
                Ready to predict
            </div>
            <div style="font-size:0.88rem;margin-top:0.4rem">
                Paste a sequence or paste a UniProt URL in the sidebar, then click
                <b>Predict Structure</b>
            </div>
        </div>
        """, unsafe_allow_html=True)

        st.markdown("#### 💡 Example proteins to try")
        examples = [
            ("Haemoglobin α-chain", "P69905", "Human oxygen carrier"),
            ("Insulin", "P01308", "Human metabolic hormone"),
            ("GFP", "P42212", "Green fluorescent protein"),
        ]
        for col, (name, acc, desc) in zip(st.columns(3), examples):
            with col:
                with st.container(border=True):
                    st.markdown(f"**{name}**")
                    st.caption(desc)
                    st.code(
                        f"https://rest.uniprot.org/uniprotkb/{acc}.fasta",
                        language=None,
                    )

    # Results
    if st.session_state.pdb_string:
        pdb_string = st.session_state.pdb_string
        plddt_df   = extract_plddt_per_residue(pdb_string)
        ss_counts  = get_secondary_structure(pdb_string)
        folded_seq = (txt or "")[:max_len]

        mean_plddt = round(float(plddt_df["plDDT"].mean()), 2) if not plddt_df.empty else None

        # Summary metrics row 
        if mean_plddt is not None:
            conf_label = (
                "Very High" if mean_plddt >= 90 else
                "Confident" if mean_plddt >= 70 else
                "Low"        if mean_plddt >= 50 else "Very Low"
            )
            c1, c2, c3, c4 = st.columns(4)
            c1.markdown(metric_card("Sequence Length", f"{len(folded_seq)} aa"), unsafe_allow_html=True)
            c2.markdown(metric_card("Mean plDDT", f"{mean_plddt}"), unsafe_allow_html=True)
            c3.markdown(metric_card("Confidence", conf_label), unsafe_allow_html=True)
            if ss_counts:
                dom  = max(ss_counts, key=ss_counts.get)
                dom_pct = round(ss_counts[dom] / sum(ss_counts.values()) * 100, 0)
                c4.markdown(
                    metric_card("Dominant SS", dom.split("/")[0], f"{dom_pct:.0f}% of residues"),
                    unsafe_allow_html=True,
                )
            st.markdown("<br>", unsafe_allow_html=True)

        # Result tabs 
        tab1, tab2, tab3, tab4 = st.tabs([
            "🔬  3D Structure",
            "📊  plDDT Analysis",
            "🧩  Secondary Structure",
            "🔤  Sequence Analysis",
        ])

        # Tab 1: 3D viewer
        with tab1:
            with st.container(border=True):
                viewer_col, info_col = st.columns([3, 1])
                with viewer_col:
                    render_mol(pdb_string, color_mode)
                with info_col:
                    st.markdown("**Color scheme**")
                    st.caption(color_mode)
                    st.markdown("---")
                    if mean_plddt is not None:
                        st.markdown("**Mean plDDT**")
                        st.markdown(plddt_badge(mean_plddt), unsafe_allow_html=True)
                    st.markdown("---")
                    st.download_button(
                        "⬇️ Download PDB",
                        data=pdb_string,
                        file_name="predicted.pdb",
                        mime="text/plain",
                        width='stretch',
                    )
                    st.caption(
                        "Open with PyMOL, ChimeraX, or "
                        "[RCSB Viewer](https://www.rcsb.org/3d-view)"
                    )

        # Tab 2: plDDT
        with tab2:
            if not plddt_df.empty:
                with st.container(border=True):
                    st.markdown(plddt_badge(mean_plddt), unsafe_allow_html=True)
                    st.markdown("<br>", unsafe_allow_html=True)

                    band_cols = st.columns(4)
                    for col, (bg, rng, label) in zip(band_cols, [
                        ("#dbeafe", "> 90",  "Very High"),
                        ("#d1fae5", "70–90", "Confident"),
                        ("#fef9c3", "50–70", "Low"),
                        ("#fee2e2", "< 50",  "Very Low"),
                    ]):
                        col.markdown(
                            f'<div style="background:{bg};border-radius:8px;padding:0.45rem;'
                            f'text-align:center;font-size:0.78rem"><b>{rng}</b><br>{label}</div>',
                            unsafe_allow_html=True,
                        )

                    st.markdown("<br>", unsafe_allow_html=True)
                    st.markdown("**Per-residue plDDT score**")
                    st.line_chart(plddt_df, width='stretch', height=300)
                    st.caption(
                        "Each point = one residue (Cα atom). "
                        "Drops below 50 often mark intrinsically disordered regions."
                    )

        # Tab 3: Secondary structure 
        with tab3:
            if ss_counts:
                with st.container(border=True):
                    total_res = sum(ss_counts.values())
                    icons = {"Alpha Helix": "🌀", "Beta Sheet": "➡️", "Coil / Loop": "〰️"}
                    m1, m2, m3 = st.columns(3)
                    for col, (name, count) in zip([m1, m2, m3], ss_counts.items()):
                        pct = round(count / total_res * 100, 1)
                        col.markdown(
                            metric_card(f"{icons[name]} {name}", str(count), f"{pct}% of residues"),
                            unsafe_allow_html=True,
                        )
                    st.markdown("<br>", unsafe_allow_html=True)
                    st.bar_chart(pd.DataFrame({"Residues": ss_counts}), width='stretch', height=280)
                    st.caption("Computed using biotite's P-SEA algorithm from the predicted PDB structure.")

        # Tab 4: Sequence analysis
        with tab4:
            if folded_seq:
                with st.container(border=True):
                    chart_col, seq_col = st.columns([2, 1])
                    with chart_col:
                        st.markdown("**Amino Acid Composition**")
                        st.bar_chart(get_aa_composition(folded_seq), width='stretch', height=320)
                    with seq_col:
                        st.markdown("**Folded sequence**")
                        wrapped = "\n".join(
                            folded_seq[i:i+10] for i in range(0, len(folded_seq), 10)
                        )
                        st.code(wrapped, language=None)
                        st.caption(
                            f"{len(folded_seq)} residues · "
                            f"{len(set(folded_seq))}/20 unique amino acids"
                        )
