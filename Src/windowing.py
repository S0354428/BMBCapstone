"""
windowing.py  --  Step 2 of the Figure 1 pipeline.

Cuts the same window out of Parent / Designed / Sequenced and prints them
stacked as a plain-text three-way alignment (DNA + protein). No figure yet.

Run:     python windowing.py
Import:  from windowing import get_window
"""

from load_data import load_all

# Mutation site, 1-based, inclusive (a607t, c608a, c609t -> codon 203)
SITE_START = 607
SITE_END = 609
FLANK = 15  # bp on each side; multiple of 3 keeps the window codon-aligned

# Row labels (placeholders, we can rename these later)
LABELS = {
    "parent": "Parent",
    "designed": "Designed",
    "sequenced": "Sequenced",
}


def get_window(data, site_start=SITE_START, site_end=SITE_END, flank=FLANK):
    """
    Return a dict with the window coordinates and, for each sequence, its DNA
    slice and its protein slice. Coordinates are 1-based.
    """
    start = site_start - flank
    end = site_end + flank
    if start < 1 or end > len(data["parent"]):
        raise ValueError("Window runs off the end of the sequence.")
    if (start - 1) % 3 != 0 or (end - start + 1) % 3 != 0:
        raise ValueError("Window is not codon-aligned; use a flank that is a multiple of 3.")

    out = {
        "start": start,
        "end": end,
        "site_start": site_start,
        "site_end": site_end,
        "first_codon": (start - 1) // 3 + 1,
        "dna": {},
        "protein": {},
    }
    for key in LABELS:
        dna = data[key][start - 1:end]
        out["dna"][key] = str(dna)
        out["protein"][key] = str(dna.translate())
    return out


def match_line(ref, other):
    """'|' where two equal-length strings agree, ' ' where they differ."""
    return "".join("|" if a == b else " " for a, b in zip(ref, other))


def protein_row(aas):
    """Spread one amino acid under the middle base of each codon."""
    return "".join(f" {aa} " for aa in aas)


def print_window(w):
    width = max(len(v) for v in LABELS.values()) + 2
    print(f"Window: nt {w['start']}-{w['end']}  "
          f"(codons {w['first_codon']}-{w['first_codon'] + len(w['protein']['parent']) - 1}), "
          f"mutation site nt {w['site_start']}-{w['site_end']}\n")

    # marker line pointing at the mutated codon
    pad = " " * width
    marks = "".join("^" if w["site_start"] <= w["start"] + i <= w["site_end"] else " "
                    for i in range(len(w["dna"]["parent"])))
    print("DNA")
    for key, label in LABELS.items():
        print(f"{label:<{width}}{w['dna'][key]}")
    print(f"{'':<{width}}{marks}\n")

    print("Match vs Parent")
    for key in ("designed", "sequenced"):
        print(f"{LABELS[key]:<{width}}{match_line(w['dna']['parent'], w['dna'][key])}")
    print()

    print("Protein")
    for key, label in LABELS.items():
        print(f"{label:<{width}}{protein_row(w['protein'][key])}")


if __name__ == "__main__":
    print_window(get_window(load_all()))
