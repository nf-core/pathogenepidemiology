process STAGE_COHORT_VCFS {
    tag "${meta.id}"
    label 'process_single'

    // Only shell symlinks are needed; reuse the concat runtime for portability.
    conda "${moduleDir}/../../nf-core/bcftools/concat/environment.yml"
    container 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'

    input:
    tuple val(meta), path(vcfs, stageAs: 'input??????.vcf.gz'), path(indices, stageAs: 'input??????.vcf.gz.tbi')

    output:
    tuple val(meta), path('ordered/*.vcf.gz', arity: '1..*'), path('ordered/*.vcf.gz.tbi', arity: '1..*'), emit: vcfs

    script:
    """
    mkdir ordered
    for input in input*.vcf.gz input*.vcf.gz.tbi; do
        ln -s ../\$input ordered/\$input
    done
    """
}
