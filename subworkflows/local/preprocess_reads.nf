
include { FASTQC       } from '../../modules/nf-core/fastqc/main'
include { BBDUK_CUSTOM as BBDUK_TRIM   } from '../../modules/local/bbduk_custom/main'
include { BBDUK_CUSTOM as BBDUK_FILTER } from '../../modules/local/bbduk_custom/main'




workflow PREPROCESS_READS {

    take:
    ch_samples
    ch_hostfasta // from workflow PREPARE_REFERENCES

    main:
    ch_adapters = channel.fromPath("${projectDir}/assets/adapters.fa", checkIfExists: true)
    ch_reads_raw = ch_samples

    // Produce QC reports per sample
    fastqc_out = FASTQC(ch_samples)

    // Trim adapters first; only the second process filters whole reads against contaminants.
    ch_reads_trimmed = BBDUK_TRIM(ch_reads_raw, ch_adapters.first())
    ch_contaminants = ch_hostfasta.first().ifEmpty { [] }
    ch_reads_clean = BBDUK_FILTER(ch_reads_trimmed.reads, ch_contaminants)

    emit:
    reads_clean   = ch_reads_clean.reads
    bbduk_stats   = ch_reads_trimmed.stats.mix(ch_reads_clean.stats)
    bbduk_logs    = ch_reads_trimmed.log.mix(ch_reads_clean.log)
    bbduk_dropped = ch_reads_clean.discarded
    trim_dropped  = ch_reads_trimmed.discarded
    samples       = ch_samples.map { meta, reads -> tuple(meta, reads, meta.instrument_platform) }
    fastqc        = fastqc_out.zip
}
