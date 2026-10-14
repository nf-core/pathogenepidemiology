
include { FASTQC       } from '../../modules/nf-core/fastqc/main'
include { BBDUK_CUSTOM } from '../../modules/local/bbduk_custom/main'




workflow PREPROCESS_READS {

    take:
    ch_samples
    ch_hostfasta // from workflow PREPARE_REFERENCES

    main:
    ch_adapters = Channel.fromPath("${projectDir}/assets/adapters.fa", checkIfExists: true)
    ch_reads_raw = ch_samples

    // Produce QC reports per sample
    fastqc_out = FASTQC(ch_samples)

    // Filter out host reads, adapters, phix
    ch_reads_clean = BBDUK_CUSTOM(ch_reads_raw, ch_hostfasta.first(), ch_adapters.first())

    emit:
    reads_clean   = ch_reads_clean.reads
    bbduk_stats   = ch_reads_clean.stats
    bbduk_logs    = ch_reads_clean.log
    bbduk_dropped = ch_reads_clean.discarded
    samples       = ch_samples.map { meta, reads -> tuple(meta, reads, meta.instrument_platform) }
    fastqc        = fastqc_out.zip
}
