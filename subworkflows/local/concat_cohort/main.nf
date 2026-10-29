include { STAGE_COHORT_VCFS } from '../../../modules/local/stage_cohort_vcfs/main'
include { BCFTOOLS_CONCAT } from '../../../modules/nf-core/bcftools/concat/main'

workflow CONCAT_COHORT {
    take:
    ch_records // [contig metadata, VCF, index], including VQSR and bypass outputs
    ch_queryfai // Reusable reference index: same contig order as its sequence dictionary.

    main:
    ch_ordered = ch_records.collect(flat: false)
        .filter { items -> !items.isEmpty() }
        .map { items ->
            def contigs = ch_queryfai.val[1].readLines().collect { row -> row.split('\t')[0] }
            def ranks = contigs.withIndex().collectEntries { contig, rank -> [(contig): rank] }
            def ids = items.collect { item -> item[0].id }
            if (ids.unique(false).size() != ids.size() || ids.any { id -> !ranks.containsKey(id) }) {
                error('Cohort VCFs must have unique contig IDs present in the query reference index')
            }
            def sorted = items.sort { left, right -> ranks[left[0].id] <=> ranks[right[0].id] }
            tuple([id: 'cohort'], sorted.collect { item -> item[1] }, sorted.collect { item -> item[2] })
        }

    // The vendored concat module sorts filenames. Numbered staged names preserve
    // reference order even when contig names and VQSR/bypass suffixes sort differently.
    STAGE_COHORT_VCFS(ch_ordered)
    BCFTOOLS_CONCAT(STAGE_COHORT_VCFS.out.vcfs)

    emit:
    vcf = BCFTOOLS_CONCAT.out.vcf
    index = BCFTOOLS_CONCAT.out.index
}
