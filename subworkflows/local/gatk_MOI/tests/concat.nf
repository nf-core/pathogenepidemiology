// Text-only integration harness: no variant analysis or biological tools.
include { PREPARE_INTERVALS } from '../../prepare_intervals/main'
include { BCFTOOLS_CONCAT } from '../../../../modules/nf-core/bcftools/concat/main'

process TEXT_VCF {
    input:
    tuple val(meta), path(bed)
    output:
    tuple val(meta), path('*.vcf.gz'), path('*.tbi')
    script:
    def prefix = "${meta.id}${meta.id == 'chr10' ? '.snp.recal' : ''}"
    """
    cp ${bed} ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}

workflow INTERVAL_CONCAT {
    take:
    bed
    fai
    main:
    PREPARE_INTERVALS(bed, channel.value(tuple([id: 'query'], file("${fai.parent}/query.fa"))),
        channel.value(tuple([id: 'query'], file("${fai.parent}/query.fa.fai"))),
        channel.value(tuple([id: 'query'], file("${fai.parent}/query.dict"))))
    TEXT_VCF(PREPARE_INTERVALS.out.intervals)
    concat_input = TEXT_VCF.out.collect(flat: false)
        .map { items ->
            tuple([id: 'cohort'], items.collect { it[1] }, items.collect { it[2] })
        }
    BCFTOOLS_CONCAT(concat_input)
    emit:
    vcf = BCFTOOLS_CONCAT.out.vcf
}
