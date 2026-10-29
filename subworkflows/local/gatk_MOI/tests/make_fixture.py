"""Deterministic, non-biological two-contig/two-sample variant-calling fixture."""

import random
from pathlib import Path

rng = random.Random(61)
references = {name: "".join(rng.choices("ACGT", k=2000)) for name in ("chr2", "chr10")}
Path("query.fa").write_text(
    "".join(f">{name}\n{seq}\n" for name, seq in references.items())
)
Path("core.bed").write_text("".join(f"{name}\t0\t2000\n" for name in references))
for name in references:
    Path(f"{name}.bed").write_text(f"{name}\t0\t2000\n")

for sample in ("sample1", "sample2"):
    header = "@HD\tVN:1.6\tSO:unsorted\n"
    header += "".join(
        f"@SQ\tSN:{name}\tLN:{len(seq)}\n" for name, seq in references.items()
    )
    header += f"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA\tLB:lib1\tPU:unit1\n"
    rows = []
    for name, reference in references.items():
        sequence = list(reference)
        sequence[599] = next(base for base in "ACGT" if base != sequence[599])
        sequence = "".join(sequence)
        for number, start in enumerate(range(300, 900, 5)):
            mate = start + 150
            qname = f"{sample}_{name}_{number}"
            for flag, position, partner, tlen in (
                (99, start, mate, 250),
                (147, mate, start, -250),
            ):
                read = sequence[position - 1 : position + 99]
                rows.append(
                    f"{qname}\t{flag}\t{name}\t{position}\t60\t100M\t=\t{partner}\t{tlen}\t{read}\t{'I' * 100}\tRG:Z:{sample}\n"
                )
    Path(f"{sample}.sam").write_text(header + "".join(rows))
