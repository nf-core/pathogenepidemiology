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
    parser.add_argument("--stub", action="store_true")
    parser.add_argument("fastas", nargs="*", type=Path)
    args = parser.parse_args()

    sample_ids = json.loads(
        base64.b64decode(args.sample_ids_base64).decode()
    )

    contigs = []
    with open(args.reference_fai) as handle:
        for line in handle:
            if line.strip():
                contigs.append(line.split("\t")[0])

    if not contigs:
        raise ValueError("Reference FAI contains no contigs")
    if len(contigs) != len(set(contigs)):
        raise ValueError("Reference FAI contains duplicate contig names")

    if not sample_ids or len(set(sample_ids)) != len(sample_ids):
        raise ValueError("Sample IDs must be present and unique")

    if args.stub:
        output_dir = Path("contigs")
        output_dir.mkdir(exist_ok=True)
        for index, contig in enumerate(contigs, start=1):
            output = output_dir / f"contig_{index:06d}.fa"
            with output.open("w") as handle:
                for sample_id in sample_ids:
                    handle.write(f">{sample_id}\nA\n")
        return

    if len(args.fastas) != len(sample_ids):
        raise ValueError(
            f"Received {len(args.fastas)} FASTA files for "
            f"{len(sample_ids)} sample IDs"
        )

    sample_records = [
        (sample_id, read_fasta(fasta))
        for sample_id, fasta in zip(sample_ids, args.fastas)
    ]

    for sample_id, records in sample_records:
        missing = set(contigs) - set(records)
        unexpected = set(records) - set(contigs)
        empty = [contig for contig in contigs if not records.get(contig)]

        if missing or unexpected or empty:
            raise ValueError(
                f"Contig mismatch for sample {sample_id}: "
                f"missing={sorted(missing)}, "
                f"unexpected={sorted(unexpected)}, "
                f"empty={sorted(empty)}"
            )

    output_dir = Path("contigs")
    output_dir.mkdir(exist_ok=True)

    for index, contig in enumerate(contigs, start=1):
        output = output_dir / f"contig_{index:06d}.fa"

        with output.open("w") as handle:
            for sample_id, records in sample_records:
                handle.write(f">{sample_id}\n")
                sequence = records[contig]

                for start in range(0, len(sequence), 80):
                    handle.write(sequence[start:start + 80] + "\n")

    print(
        f"Created {len(contigs)} contig FASTAs "
        f"for {len(sample_ids)} samples"
    )


if __name__ == "__main__":
    main()
