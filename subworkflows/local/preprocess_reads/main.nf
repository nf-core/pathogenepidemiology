include { FASTQC } from '../../../modules/nf-core/fastqc/main'
include { BBDUK_CUSTOM as BBDUK_TRIM } from '../../../modules/local/bbduk_custom/main'
include { BBDUK_CUSTOM as BBDUK_FILTER } from '../../../modules/local/bbduk_custom/main'
include { NANOPLOT } from '../../../modules/nf-core/nanoplot/main'
include { FILTLONG } from '../../../modules/nf-core/filtlong/main'




workflow PREPROCESS_READS {
    take:
    ch_samples
    ch_hostfasta // reusable path value, or [] to filter only against PhiX

    main:
    def ch_by_platform = ch_samples.branch { meta, _reads ->
        ont: meta.instrument_platform == 'OXFORD_NANOPORE'
        other: true
    }

    ch_adapters = channel.value(file("${projectDir}/assets/adapters.fa", checkIfExists: true))
    fastqc_out = FASTQC(ch_by_platform.other)
    ch_reads_trimmed = BBDUK_TRIM(ch_by_platform.other, ch_adapters)
    ch_reads_clean = BBDUK_FILTER(ch_reads_trimmed.reads, ch_hostfasta)

    def ch_raw_ont = ch_by_platform.ont.map { meta, reads ->
        if (!meta.single_end || reads.size() != 1) {
            error("Oxford Nanopore sample ${meta.id} must provide one FASTQ in fastq_1 and leave fastq_2 empty.")
        }
        tuple(meta, reads[0])
    }
    nanoplot_out = NANOPLOT(ch_raw_ont)
    filtlong_out = FILTLONG(ch_raw_ont.map { meta, read -> tuple(meta, [], read) })

    emit:
    reads_clean = ch_reads_clean.reads
    bbduk_stats = ch_reads_trimmed.stats.mix(ch_reads_clean.stats)
    bbduk_logs = ch_reads_trimmed.log.mix(ch_reads_clean.log)
    bbduk_dropped = ch_reads_clean.discarded
    fastqc = fastqc_out.zip
    ont_reads = filtlong_out.reads
    nanoplot = nanoplot_out.txt
    filtlong_log = filtlong_out.log
}
