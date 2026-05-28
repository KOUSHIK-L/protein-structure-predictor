# ESMFold Protein Structure Prediction

A Streamlit web app that predicts 3D protein structure from an amino acid sequence using **Meta AI's ESMFold** — no alignment required, results in seconds.

## Live Demo

[![Streamlit App](https://img.shields.io/badge/Streamlit-FF4B4B?style=for-the-badge&logo=Streamlit&logoColor=white)](https://protein-structure-predictor-using-esm-fold.streamlit.app/)

---

## Features

- **Three input modes** — paste a raw sequence, provide a UniProt FASTA URL, or upload a `.fasta` / `.txt` file
- **UniProt metadata** — when a UniProt URL is provided, the app automatically fetches protein name, gene, organism, and functional annotation
- **Interactive 3D viewer** — rotate and zoom the predicted structure with three color schemes:
  - Spectrum (rainbow by residue position)
  - plDDT Confidence (blue = high, red = low — AlphaFold convention)
  - Secondary Structure (Jmol coloring)
- **Per-residue plDDT chart** — line chart of confidence scores along the sequence, highlighting disordered regions
- **Secondary structure breakdown** — alpha helix, beta sheet, and coil/loop residue counts computed via biotite's P-SEA algorithm
- **Amino acid composition** — bar chart of all 20 standard amino acids in the input sequence
- **Sequence length guard** — warns when sequences exceed 400 aa (ESMFold's reliable range)
- **Result caching** — repeated predictions for the same sequence skip the API call
- **PDB download** — save the predicted structure for use in PyMOL, ChimeraX, or other tools

---

## Understanding the Output

| plDDT | Interpretation |
|-------|----------------|
| > 90  | Very high confidence — backbone is reliable |
| 70–90 | Confident — generally accurate |
| 50–70 | Low confidence — treat with caution |
| < 50  | Likely disordered or unreliable |

---

## Running Locally

```bash
# 1. Clone the repo
git clone https://github.com/KOUSHIK-L/protein-structure-predictor.git
cd protein-structure-predictor

# 2. Create and activate a virtual environment
python -m venv venv
source venv/bin/activate      # Windows: venv\Scripts\activate

# 3. Install dependencies
pip install -r requirements.txt

# 4. Run the app
streamlit run streamlit_app.py
```

---

## Dependencies

| Package | Purpose |
|---------|---------|
| `streamlit` | Web app framework |
| `streamlit-option-menu` | Navigation menu |
| `stmol` | Streamlit + py3Dmol bridge |
| `py3Dmol` | 3D molecular viewer |
| `biotite` | Secondary structure annotation (P-SEA) |
| `numpy` | Numerical operations |
| `pandas` | Data tables and charts |
| `requests` | ESMFold API + UniProt API calls |

---

## How It Works

1. User provides a protein sequence (raw, FASTA URL, or file)
2. The app sends the sequence to the [ESMFold REST API](https://esmatlas.com/about)
3. The API returns a PDB-format structure predicted by ESM-2 + folding head
4. The app parses B-factor fields for per-residue plDDT scores
5. Biotite's `annotate_sse` (P-SEA algorithm) classifies residues into helix, sheet, or coil
6. Everything is visualised in-browser with py3Dmol and Streamlit charts

---

## Limitations

- ESMFold works best on sequences ≤ 400 amino acids
- The public ESMFold API (`api.esmatlas.com`) can be intermittently unavailable; as an alternative, use [ColabFold](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)
- ESMFold does not use multiple sequence alignments, so accuracy on highly divergent or intrinsically disordered proteins may be lower than AlphaFold2

---

## References

- Lin et al. (2023). *Evolutionary-scale prediction of atomic-level protein structure with a language model.* Science. [DOI: 10.1126/science.ade2574](https://www.science.org/doi/10.1126/science.ade2574)
- [ESM Metagenomic Atlas — Meta AI blog](https://ai.facebook.com/blog/protein-folding-esmfold-metagenomics/)

---

## Author

**Koushik Loganathan**  
M.Sc. Life Science Informatics, University of Bonn  
[GitHub](https://github.com/KOUSHIK-L) · [Live App](https://protein-structure-predictor-using-esm-fold.streamlit.app/)
