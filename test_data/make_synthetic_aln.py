"""Build a tiny synthetic MSA from a reference GenBank for smoke-testing.

We replicate the reference sequence N times and sprinkle random SNPs, with
a couple of higher-variability "hotspot" windows to exercise the ranking
logic. Output: test_data/synthetic.fasta.
"""

from __future__ import annotations

import random
from pathlib import Path

from Bio import SeqIO

HERE = Path(__file__).parent
REF_GB = HERE / "Salvia_chinensis.gb"
OUT = HERE / "synthetic.fasta"

N_COPIES = 6
BASE_SNP_RATE = 0.004
HOTSPOT_WINDOWS = [(40_000, 41_200), (95_000, 96_500), (130_000, 131_000)]
HOTSPOT_RATE = 0.06
RNG = random.Random(42)


def mutate(seq: str, rng: random.Random) -> str:
    bases = list(seq.upper())
    alphabets = {"A": "CGT", "C": "AGT", "G": "ACT", "T": "ACG"}
    L = len(bases)

    hotspot_mask = [False] * L
    for s, e in HOTSPOT_WINDOWS:
        for i in range(max(0, s), min(L, e)):
            hotspot_mask[i] = True

    for i, b in enumerate(bases):
        if b not in alphabets:
            continue
        rate = HOTSPOT_RATE if hotspot_mask[i] else BASE_SNP_RATE
        if rng.random() < rate:
            bases[i] = rng.choice(alphabets[b])
    return "".join(bases)


def main():
    rec = next(SeqIO.parse(str(REF_GB), "genbank"))
    ref = str(rec.seq).upper()
    with OUT.open("w") as fh:
        fh.write(f">{rec.id}\n{ref}\n")
        for i in range(N_COPIES):
            rng = random.Random(100 + i)
            mut = mutate(ref, rng)
            fh.write(f">sim_{i+1}\n{mut}\n")
    print(f"Wrote {OUT}  ({N_COPIES + 1} sequences, length {len(ref):,})")


if __name__ == "__main__":
    main()
