process SPLIT_CONSENSUS_CONTIGS {

    tag 'split_consensus_contigs'
    label 'process_single'

    container 'docker.io/library/python:3.12-bookworm'

    input:
    tuple val(sample_ids), path(fastas, stageAs: 'sample???.fa'), path(reference_fai)

    output:
    path 'contigs/*.fa', arity: '1..*', emit: contigs

    script:
    def sample_ids_json = groovy.json.JsonOutput.toJson(sample_ids)
    def sample_ids_b64 = sample_ids_json.bytes.encodeBase64().toString()

    """
    python3 - <<'PY'
    import base64
    import json
    from pathlib import Path

    sample_ids = json.loads(base64.b64decode("${sample_ids_b64}"))
    fasta_files = [
        Path(f"sample{i + 1:03d}.fa")
        for i in range(len(sample_ids))
    ]

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

    contigs = []
    with open("${reference_fai}") as handle:
        for line in handle:
            if line.strip():
                contigs.append(line.split("\\t")[0])

    if not contigs:
        raise ValueError("Reference FAI contains no contigs")
    if len(contigs) != len(set(contigs)):
        raise ValueError("Reference FAI contains duplicate contig names")

    if not sample_ids or len(set(sample_ids)) != len(sample_ids):
        raise ValueError("Sample IDs must be present and unique")

    if len(fasta_files) != len(sample_ids):
        raise ValueError(
            f"Received {len(fasta_files)} FASTA files for "
            f"{len(sample_ids)} sample IDs"
        )

    sample_records = [
        (sample_id, read_fasta(fasta))
        for sample_id, fasta in zip(sample_ids, fasta_files)
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
                handle.write(f">{sample_id}\\n")
                sequence = records[contig]

                for start in range(0, len(sequence), 80):
                    handle.write(sequence[start:start + 80] + "\\n")

    print(
        f"Created {len(contigs)} contig FASTAs "
        f"for {len(sample_ids)} samples"
    )
    PY
    """
    stub:
    def sample_ids_json = groovy.json.JsonOutput.toJson(sample_ids)
    def sample_ids_b64 = sample_ids_json.bytes.encodeBase64().toString()

    """
    mkdir -p contigs

    python3 - <<'PY'
    import base64
    import json
    from pathlib import Path

    sample_ids = json.loads(
        base64.b64decode("${sample_ids_b64}").decode()
    )

    with open("${reference_fai}") as fai:
        contigs = [line.split("\\t")[0] for line in fai if line.strip()]

    for index, contig in enumerate(contigs, start=1):
        output = Path(f"contigs/contig_{index:06d}.fa")
        with output.open("w") as handle:
            for sample_id in sample_ids:
                handle.write(f">{sample_id}\\nA\\n")
    PY
    """
}
