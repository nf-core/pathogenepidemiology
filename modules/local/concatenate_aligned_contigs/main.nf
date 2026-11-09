process CONCATENATE_ALIGNED_CONTIGS {

    tag 'concatenate_aligned_contigs'
    label 'process_single'

    container 'docker.io/library/python:3.12-bookworm'

    input:
    tuple val(sample_ids), path(alignments), path(reference_fai)

    output:
    path 'cohort.concatenated.fas', emit: alignment

    script:
    def sample_ids_json = groovy.json.JsonOutput.toJson(sample_ids)
    def sample_ids_b64 = sample_ids_json.bytes.encodeBase64().toString()

    """
    python3 - <<'PY'
    import base64
    import json
    from pathlib import Path

    sample_ids = json.loads(base64.b64decode("${sample_ids_b64}"))

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

    if not sample_ids or len(set(sample_ids)) != len(sample_ids):
        raise ValueError("Sample IDs must be present and unique")

    contigs = []
    with open("${reference_fai}") as handle:
        for line in handle:
            if line.strip():
                contigs.append(line.split("\\t")[0])

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
            output.write(f">{sample_id}\\n")
            sequence = "".join(concatenated[sample_id])
            for start in range(0, len(sequence), 80):
                output.write(sequence[start:start + 80] + "\\n")

    print(
        f"Concatenated {len(contigs)} contigs "
        f"for {len(sample_ids)} samples"
    )
    PY
    """

    stub:
    """
    touch cohort.concatenated.fas
    """
}
