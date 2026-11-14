# BMBCapstone

GFP variant design and sequencing analysis pipeline.

- `Src/load_data.py` – load parent, designed, and sequenced GFP; sanity checks
- `Src/make_figure.py` – generate Figure 1
- `Data/` – input sequences

## Setup

    python3 -m venv .venv
    source .venv/bin/activate
    pip install -r Src/requirements.txt

## Run

    cd Src
    python load_data.py
    python make_figure.py
