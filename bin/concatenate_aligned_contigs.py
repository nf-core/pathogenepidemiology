#!/usr/bin/env python3

import argparse
import base64
import json
from pathlib import Path


def read_fasta(path):
    records = {}
    current = None
    sequences = []

    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    records[current] = "".join(sequences)
                current = line[1:].split()[0]
                if current in records:
                    raise ValueError(
                        f"Duplicate FASTA record {current} in {path}"
                    )
                sequences = []
            else:
                if current is None:
                    raise ValueError(
                        f"Sequence before FASTA header in {path}"
                    )
                sequences.append(line)

    if current is not None:
        records[current] = "".join(sequences)

    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-ids-base64", required=True)
    parser.add_argument("--reference-fai", required=True)
    args = parser.parse_args()

    sample_ids = json.loads(
        base64.b64decode(args.sample_ids_base64).decode()
    )

    if not sample_ids or len(set(sample_ids)) != len(sample_ids):
        raise ValueError("Sample IDs must be present and unique")

    contigs = []
    with open(args.reference_fai) as handle:
        for line in handle:
            if line.strip():
                contigs.append(line.split("\t")[0])

    if not contigs:
        raise ValueError("Reference FAI contains no contigs")
    if len(contigs) != len(set(contigs)):
        raise ValueError("Reference FAI contains duplicate contig names")

    expected_files = {
        f"contig_{index:06d}.fas"
        for index in range(1, len(contigs) + 1)
    }
    observed_files = {
        path.name for path in Path(".").glob("contig_*.fas")
    }

    if observed_files != expected_files:
        raise ValueError(
            f"Alignment file mismatch: "
            f"missing={sorted(expected_files - observed_files)}, "
            f"unexpected={sorted(observed_files - expected_files)}"
        )

    concatenated = {sample_id: [] for sample_id in sample_ids}

    for index, contig in enumerate(contigs, start=1):
        alignment = Path(f"contig_{index:06d}.fas")
        records = read_fasta(alignment)

        missing = set(sample_ids) - set(records)
        unexpected = set(records) - set(sample_ids)
        empty = [sample for sample in sample_ids if not records.get(sample)]

        if missing or unexpected or empty:
            raise ValueError(
                f"Sample mismatch in {alignment} for contig {contig}: "
                f"missing={sorted(missing)}, "
                f"unexpected={sorted(unexpected)}, "
                f"empty={sorted(empty)}"
            )

        lengths = {len(sequence) for sequence in records.values()}
        if len(lengths) != 1:
            raise ValueError(
                f"Sequences in {alignment} have unequal alignment lengths"
            )

        for sample_id in sample_ids:
            concatenated[sample_id].append(records[sample_id])

    with open("cohort.concatenated.fas", "w") as output:
        for sample_id in sample_ids:
            output.write(f">{sample_id}\n")
            sequence = "".join(concatenated[sample_id])
            for start in range(0, len(sequence), 80):
                output.write(sequence[start:start + 80] + "\n")

    print(
        f"Concatenated {len(contigs)} contigs "
        f"for {len(sample_ids)} samples"
    )


if __name__ == "__main__":
    main()
