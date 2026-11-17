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
    concatenate_aligned_contigs.py \\
        --sample-ids-base64 '${sample_ids_b64}' \\
        --reference-fai '${reference_fai}'
    """

    stub:
    """
    touch cohort.concatenated.fas
    """
}
