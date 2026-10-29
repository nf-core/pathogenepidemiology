include { GATK_MOI } from '../main.nf'
include { CHECK_COHORT } from '../../concat_cohort/tests/main.nf'

process MAKE_COHORT_FIXTURE {
    label 'process_single'
    container 'community.wave.seqera.io/library/gatk4_gcnvkernel_htslib_samtools:d3becb6465454c35'

    input:
    path generator

    output:
    path 'query.fa', emit: fasta
    path 'query.fa.fai', emit: fai
    path 'sample*.bam', emit: bams
    path 'core.bed', emit: core
    path 'chr2.bed', emit: chr2
    path 'chr10.bed', emit: chr10

    script:
    """
    python ${generator}
    samtools faidx query.fa
    for sample in sample1 sample2; do
        samtools sort -o \$sample.bam \$sample.sam
    done
    """
}

workflow REAL_GATK_COHORT {
    main:
    MAKE_COHORT_FIXTURE(file("${projectDir}/subworkflows/local/gatk_MOI/tests/make_fixture.py"))
    ch_bams = MAKE_COHORT_FIXTURE.out.bams.flatMap { bams -> bams.collect { bam -> tuple([id: bam.baseName], bam) } }
    ch_fasta = MAKE_COHORT_FIXTURE.out.fasta.map { fasta -> tuple([id: 'query'], fasta) }.first()
    ch_fai = MAKE_COHORT_FIXTURE.out.fai.map { fai -> tuple([id: 'query'], fai) }.first()
    ch_intervals = MAKE_COHORT_FIXTURE.out.chr10.map { bed -> tuple([id: 'chr10'], bed) }
        .mix(MAKE_COHORT_FIXTURE.out.chr2.map { bed -> tuple([id: 'chr2'], bed) })
    GATK_MOI(ch_bams, ch_fasta, ch_fai, channel.value([]), channel.value([]), MAKE_COHORT_FIXTURE.out.core.first(), ch_intervals)
    CHECK_COHORT(GATK_MOI.out.varcalls_s)

    emit:
    records = CHECK_COHORT.out.records
    samples = CHECK_COHORT.out.samples
    indexed = CHECK_COHORT.out.indexed
    cohort = GATK_MOI.out.varcalls_s
    qc = GATK_MOI.out.insert_size_metrics
}
