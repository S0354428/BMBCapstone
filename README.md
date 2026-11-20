# BMBCapstone

GFP variant design, sequencing analysis, and protein quantification pipeline.

## Figure 1: Sequencing comparison

- `Src/load_data.py` – load parent, designed, and sequenced GFP; sanity checks
- `Src/windowing.py` – cut the same window around the mutation site from all three sequences
- `Src/make_figure.py` – generate Figure 1
- `Data/` inputs: `parent.md`, `designed.md`, and the sequencing read (`11201415235-1_1_1-1.fasta`)

## Figure 2: Bradford protein quantification

- `Src/load_bradford.py` – load BSA standards and Plate 2, map wells to samples
- `Src/fit_standards.py` – blank-correct the standards and fit the linear standard curve
- `Src/quantify.py` – convert absorbance to mg/mL, apply dilution factors, write `Results/` CSVs
- `Src/make_figure2.py` – generate Figure 2
- `Data/` inputs: `bradford_exploratory_selected_triplicates.csv` (standards) and `Bradford__Plate_2_-_Sheet1_.csv` (samples)

## Outputs

- `Figures/` – `figure1` and `figure2` as PNG, PDF, and SVG
- `Results/` – protein concentrations for every sample and dilution, plus the reported value per sample

## Setup

    python3 -m venv .venv
    source .venv/bin/activate
    pip install -r Src/requirements.txt

## Run

Run everything from the `Src` folder.

Figure 1:

    cd Src
    python load_data.py
    python make_figure.py

Figure 2 (each of the first three scripts prints a check you can read):

    cd Src
    python load_bradford.py
    python fit_standards.py
    python quantify.py
    python make_figure2.py
