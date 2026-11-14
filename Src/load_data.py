"""
load_data.py  --  Step 1 of the Figure 1 pipeline.

Loads the three sequences, puts them on the same strand and the same
coordinates, and prints a sanity check. No figure code lives here.

    Parent     : parental GFP template (what we started with)
    Designed   : expected mutant (parent + the primer's changes)
    Sequenced  : the GFP region pulled out of the full-plasmid sequencing read

Run from anywhere:   python src/load_data.py
Import elsewhere:    from load_data import load_all
"""

from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

# --------------------------------------------------------------------------
# Paths  (src/ and Data/ are siblings inside BMBCap/)
# --------------------------------------------------------------------------
DATA_DIR = Path(__file__).resolve().parent.parent / "Data"

PARENT_FILE = DATA_DIR / "parent.md"
DESIGNED_FILE = DATA_DIR / "designed.md"
SEQUENCED_FILE = DATA_DIR / "11201415235-1_1_1-1.fasta"

# Mutagenic primer pair, 5' -> 3' (used only as a cross-check)
PRIMER_FWD = "CGACAACCACTACCTGAGCTATCAGTCCAAGCTGAGCAAAG"
PRIMER_REV = "CTTTGCTCAGCTTGGACTGATAGCTCAGGTAGTGGTTGTCG"


# --------------------------------------------------------------------------
# Reading
# --------------------------------------------------------------------------
def read_seq(path):
    """Read a one-record FASTA file (any extension) as an uppercase Seq."""
    record = SeqIO.read(path, "fasta")
    return Seq(str(record.seq).upper().replace(" ", "").replace("\n", ""))


# --------------------------------------------------------------------------
# Pulling GFP out of the plasmid read
# --------------------------------------------------------------------------
def extract_gfp(read, parent):
    """
    Find the parent GFP inside a full-plasmid read and return it as a Seq of
    the same length, on the parent's strand.

    - Tries both strands of the read.
    - Searches a doubled copy of the read, because the plasmid is circular
      and the GFP can wrap across the end of the file.
    - Uses a local aligner (not an exact string search) so the GFP is still
      found if the read carries a mutation or a sequencing error.

    Returns (gfp_seq, info_dict).
    """
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -3
    aligner.open_gap_score = -6
    aligner.extend_gap_score = -2

    best = None
    for strand, candidate in (("forward", read), ("reverse complement", read.reverse_complement())):
        doubled = candidate + candidate  # handle wrap-around
        aln = aligner.align(doubled, parent)[0]
        if best is None or aln.score > best["score"]:
            best = {"score": aln.score, "strand": strand, "seq": doubled, "aln": aln}

    t_blocks, q_blocks = best["aln"].aligned
    # Where parent base 0 would sit in the doubled read, even if the local
    # alignment stopped a few bases short of the edge.
    start = int(t_blocks[0][0] - q_blocks[0][0])
    if start < 0 or start + len(parent) > len(best["seq"]):
        raise ValueError("GFP not fully contained in the read; check the input files.")

    gfp = best["seq"][start:start + len(parent)]
    info = {
        "strand": best["strand"],
        "start_in_read": start % len(read),  # 0-based, on the strand shown
        "wraps_around_end": (start % len(read)) + len(parent) > len(read),
        "score": best["score"],
    }
    return gfp, info


# --------------------------------------------------------------------------
# Comparison helpers
# --------------------------------------------------------------------------
def mismatches(a, b):
    """0-based positions where two equal-length sequences differ."""
    if len(a) != len(b):
        raise ValueError(f"Length mismatch: {len(a)} vs {len(b)}")
    return [i for i, (x, y) in enumerate(zip(a, b)) if x != y]


def protein_changes(a, b):
    """(1-based residue, aa_in_a, aa_in_b) wherever the translations differ."""
    usable = len(a) - (len(a) % 3)  # trim trailing partial codon
    pa = str(a[:usable].translate())
    pb = str(b[:usable].translate())
    return [(i + 1, x, y) for i, (x, y) in enumerate(zip(pa, pb)) if x != y]


# --------------------------------------------------------------------------
# Main loader
# --------------------------------------------------------------------------
def load_all():
    """Return {'parent': Seq, 'designed': Seq, 'sequenced': Seq, 'read_info': dict}."""
    parent = read_seq(PARENT_FILE)
    designed = read_seq(DESIGNED_FILE)
    read = read_seq(SEQUENCED_FILE)
    sequenced, read_info = extract_gfp(read, parent)
    read_info["full_read_length"] = len(read)
    return {"parent": parent, "designed": designed, "sequenced": sequenced, "read_info": read_info}


# --------------------------------------------------------------------------
# Sanity check
# --------------------------------------------------------------------------
def sanity_check(data):
    parent, designed, sequenced = data["parent"], data["designed"], data["sequenced"]
    info = data["read_info"]

    print("=== Loaded sequences ===")
    print(f"Parent    : {len(parent)} nt")
    print(f"Designed  : {len(designed)} nt")
    print(f"Sequenced : {len(sequenced)} nt  (from a {info['full_read_length']} bp read)")
    print(f"  GFP found on the {info['strand']} strand, starting at read position "
          f"{info['start_in_read'] + 1} (1-based); wraps past end of file: {info['wraps_around_end']}")

    print("\n=== Parent vs Designed ===")
    d_pos = mismatches(parent, designed)
    print("DNA differences (1-based):",
          [(p + 1, parent[p], designed[p]) for p in d_pos] or "none")
    print("Protein changes:", protein_changes(parent, designed) or "none")

    print("\n=== Parent vs Sequenced ===")
    s_pos = mismatches(parent, sequenced)
    print("DNA differences (1-based):",
          [(p + 1, parent[p], sequenced[p]) for p in s_pos] or "none")
    print("Protein changes:", protein_changes(parent, sequenced) or "none")

    print("\n=== Designed vs Sequenced ===")
    ds_pos = mismatches(designed, sequenced)
    print("DNA differences (1-based):",
          [(p + 1, designed[p], sequenced[p]) for p in ds_pos] or "none")
    print("Mutation incorporated?", "YES" if not ds_pos else "NO")

    print("\n=== Primer cross-check ===")
    print("Forward primer found in Designed:", PRIMER_FWD in str(designed))
    print("Reverse primer is reverse complement of forward:",
          str(Seq(PRIMER_REV).reverse_complement()) == PRIMER_FWD)
    print("Forward primer found in Sequenced:", PRIMER_FWD in str(sequenced))


if __name__ == "__main__":
    sanity_check(load_all())
