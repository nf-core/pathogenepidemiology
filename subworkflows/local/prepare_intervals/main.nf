include { GATK4_SPLITINTERVALS } from '../../../modules/nf-core/gatk4/splitintervals/main'
include { GATK4_INTERVALLISTTOBED } from '../../../modules/nf-core/gatk4/intervallisttobed/main'

def validateLocalBed(value) {
    if (!(value instanceof CharSequence) || !(value ==~ /.+\.bed/) || value.contains('://')) {
        error '--intervals must be a single local uncompressed BED file (.bed).'
    }
    def bed = file(value, checkIfExists: true, glob: false)
    if (!bed.isFile()) {
        error "--intervals must be a single file: ${value}"
    }
    if (bed.toFile().length() == 0) {
        error "--intervals BED file must not be empty: ${value}"
    }
    return bed
}

workflow PREPARE_INTERVALS {
    take:
    intervals
    ch_queryfasta
    ch_queryfai
    ch_querydict

    main:
    // Select once: the original BED also drives core-read filtering and coverage.
    def bed = validateLocalBed(intervals != null
        ? intervals
        : "${projectDir}/assets/Pf3D7_core.bed")
    // The upstream module interpolates staged basenames into its command.
    if (!(bed.name ==~ /[A-Za-z0-9_][A-Za-z0-9._-]*/)) {
        error '--intervals BED basename must contain only letters, digits, dots, underscores or hyphens and start with a letter, digit or underscore'
    }
    core_bed = channel.value(bed)
    GATK4_SPLITINTERVALS(core_bed.map { tuple([id: 'intervals'], it) },
        ch_queryfasta, ch_queryfai, ch_querydict)
    ch_interval_lists = GATK4_SPLITINTERVALS.out.split_intervals.flatMap { meta, interval_files ->
        // Nextflow returns a Path rather than a list for a single match.
        def files = interval_files instanceof List ? interval_files : [interval_files]
        files.collect { interval ->
            def rows = interval.readLines().findAll { it && !it.startsWith('@') }
            def contigs = rows.collect { it.split('\t')[0] }.unique()
            if (contigs.size() != 1) {
                error "Expected one contig in ${interval}, found ${contigs}"
            }
            def contig = contigs.first()
            if (!(contig ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) {
                error "Unsupported contig '${contig}' for output filenames in ${interval}"
            }
            tuple([id: contig], interval)
        }
    }
        .collect(flat: false)
        .flatMap { interval_outputs ->
            // Validate all outputs before conversion or downstream sample scatter.
            interval_outputs.groupBy { it[0].id }.each { contig, outputs ->
                if (outputs.size() > 1) {
                    error "Expected one interval file per contig; duplicate contig '${contig}' in ${outputs.collect { it[1].name }}"
                }
            }
            interval_outputs
        }
    GATK4_INTERVALLISTTOBED(ch_interval_lists)
    ch_intervals = GATK4_INTERVALLISTTOBED.out.bed

    emit:
    core_bed
    intervals = ch_intervals
}
