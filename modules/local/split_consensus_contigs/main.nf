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
    split_consensus_contigs.py \\
        --sample-ids-base64 '${sample_ids_b64}' \\
        --reference-fai '${reference_fai}' \\
        ${fastas.join(' ')}
    """

    stub:
    def sample_ids_json = groovy.json.JsonOutput.toJson(sample_ids)
    def sample_ids_b64 = sample_ids_json.bytes.encodeBase64().toString()

    """
    split_consensus_contigs.py \\
        --sample-ids-base64 '${sample_ids_b64}' \\
        --reference-fai '${reference_fai}' \\
        --stub
    """
}
