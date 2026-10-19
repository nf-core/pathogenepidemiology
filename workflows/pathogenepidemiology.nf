/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_REFERENCES                   } from '../subworkflows/local/prepare_references'
include { PREPROCESS_READS                     } from '../subworkflows/local/preprocess_reads'
include { MINIMAP2_INDEX                       } from '../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN                       } from '../modules/nf-core/minimap2/align/main'
include { BWAMEM3_INDEX                        } from '../modules/nf-core/bwamem3/index/main'
include { BWAMEM3_MEM                          } from '../modules/nf-core/bwamem3/mem/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_MM2 } from '../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_BM3 } from '../modules/nf-core/samtools/stats/main'
include { GATK_MOI                             } from '../subworkflows/local/gatk_MOI/main'
include { PHYLOGENY                            } from '../subworkflows/local/phylogeny/phylogeny.nf'
include { BCFTOOLS_STATS                       } from '../modules/nf-core/bcftools/stats/main'
//include { HMMIBDRS                             } from '../modules/local/hmmibdrs/main'
include { MULTIQC                              } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                     } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                 } from '../subworkflows/nf-core/utils_nfcore_pipeline/main'
include { softwareVersionsToYAML               } from '../subworkflows/nf-core/utils_nfcore_pipeline/main'
include { methodsDescriptionText               } from '../subworkflows/local/utils_nfcore_pathogenepidemiology_pipeline/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PATHOGENEPIDEMIOLOGY {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // SUBWORKFLOW: Download and index reference fasta files
    //
    PREPARE_REFERENCES(
        params.queryurl,
        params.hosturl
    )
    def ch_queryfasta = PREPARE_REFERENCES.out.queryfasta
    def ch_hostfasta  = PREPARE_REFERENCES.out.hostfasta
    def ch_queryfai   = PREPARE_REFERENCES.out.queryfai

    //
    // SUBWORKFLOW: FastQC on reads, BBDUK adapter/contaminant removal
    //
    PREPROCESS_READS(
        ch_samplesheet,
        ch_hostfasta,
    )
    def ch_reads_clean   = PREPROCESS_READS.out.reads_clean
    def ch_bbduk_stats   = PREPROCESS_READS.out.bbduk_stats
    def ch_bbduk_logs    = PREPROCESS_READS.out.bbduk_logs
    def ch_bbduk_dropped = PREPROCESS_READS.out.bbduk_dropped
    def ch_samples       = PREPROCESS_READS.out.samples

    ch_multiqc_files = ch_multiqc_files.mix(
        PREPROCESS_READS.out.fastqc.map { _meta, files -> files },
        ch_bbduk_stats.map { _meta, files -> files },
        ch_bbduk_logs.map { _meta, files -> files },
        ch_bbduk_dropped.map { _meta, files -> files }
    )

    //
    // LONG READ TRACK - ALIGNMENT + VAR CALLING
    //

    // Filter short reads out
    def ch_samples_l = ch_reads_clean
        .join(ch_samples.map { meta, _reads, platform -> [ meta, platform ] }, by: 0)
        .filter { _meta, _reads, platform -> platform == 'OXFORD_NANOPORE' }
        .map { meta, reads, _platform -> [ meta, reads ] }

    // Alignment of long reads
    def minimap2_index = MINIMAP2_INDEX(ch_queryfasta)
    def aligned_l = MINIMAP2_ALIGN(
        ch_samples_l,
        minimap2_index.index.first(),   // reference as tuple
        true,                           // bam_format
        'bai',                          // bam_index_extension
        false,                          // cigar_paf_format
        false                           // cigar_bam
    )

    // Minimap2 stats input
    def ch_mm2stats_input = aligned_l.bam
        .join(aligned_l.index, by: 0)
        .map { meta, bam, bai -> [ meta, bam, bai ] }

    // Variant calling of long reads
    //varcalls_l = // CLAIR3 removed, to be replaced as per issue #4

    //
    // SHORT READ TRACK - ALIGNMENT + VAR CALLING
    //

    // Filter long reads out
    def ch_samples_s = ch_reads_clean
        .join(ch_samples.map { meta, _reads, platform -> [ meta, platform ] }, by: 0)
        .filter { _meta, _reads, platform -> platform == 'ILLUMINA' || platform == 'DNBSEQ' }
        .map { meta, reads, _platform -> [ meta, reads ] }

    // Alignment of short reads
    def bwamem3_index = BWAMEM3_INDEX(ch_queryfasta)
    def aligned_s = BWAMEM3_MEM(
        ch_samples_s,
        bwamem3_index.index.first(),
        ch_queryfasta.first(),
        true                            // sort the bam for gatk
    )

    // BWA-MEM3 stats input
    def ch_bm3stats_input = aligned_s.aligned
        .join(aligned_s.index, by: 0)
        .map { meta, bam, bai -> [ meta, bam, bai ] }

    def ch_wgs_bam_s      = aligned_s.aligned.filter { meta, _bam -> meta.library_strategy == 'WGS' }
    // Amplicon BAMs — reserved for the bcftools-based caller (not yet wired)
    def ch_amplicon_bam_s = aligned_s.aligned.filter { meta, _bam -> meta.library_strategy != 'WGS' }

    def gatk_out = GATK_MOI(
        ch_wgs_bam_s,
        ch_queryfasta,
        ch_queryfai,
        Channel.fromPath("${projectDir}/assets/Strains.2kb.vcf.gz", checkIfExists: true),
        Channel.fromPath("${projectDir}/assets/Strains.2kb.vcf.gz.tbi", checkIfExists: true)
    )
    def varcalls_s = gatk_out.varcalls_s

    // Generate one consensus FASTA per WGS sample from the cohort BCF
    def ch_reference = ch_queryfasta.map { _meta, fasta -> fasta }.first()
    def ch_consensus_input = varcalls_s
        .combine(ch_wgs_bam_s.map { meta, _bam -> meta })
        .combine(ch_reference)
        .map { cohort_meta, bcf, csi, sample_meta, fasta ->
            [ sample_meta, bcf, csi, fasta, [] ]
        }
    def phylogeny_s = PHYLOGENY(ch_consensus_input)

    def ch_bcftools_stats = BCFTOOLS_STATS(
        varcalls_s,        // [meta, bcf, csi]
        [[:], []],         // regions
        [[:], []],         // targets
        [[:], []],         // samples
        [[:], []],         // exons
        [[:], []]          // fasta
    )


    // TODO: MERGE varcalls_l AND varcalls_s CHANNELS // maybe fixed because a new channel can be made from a dir with all
/*
    varcalls = varcalls_s

    // Filter to sites where every sample has a populated AD.
    // GATK emits "." at low-coverage or filtered positions, which hmmibd-rs
    // rejects at parse time with NumericaValueEmptyInt.
    ch_hmmibd_input = BCFTOOLS_VIEW(
        varcalls,
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)

    hmmibd_out = HMMIBDRS(
        ch_hmmibd_input.vcf.join(ch_hmmibd_input.index, by: 0)
    )
*/

    //
    // POST PROCESSING
    //

    // samtools stats reference: tuple val(meta2), path(fasta), path(fai)
    def statsrefs = ch_queryfasta
        .join(ch_queryfai, by: 0)
        .map { meta, ref, fai -> [ meta, ref, fai ] }

    def ch_mm2stats = SAMTOOLS_STATS_MM2(ch_mm2stats_input, statsrefs.first())
    def ch_bm3stats = SAMTOOLS_STATS_BM3(ch_bm3stats_input, statsrefs.first())

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_mm2stats.stats.map { _meta, files -> files },
        ch_bm3stats.stats.map { _meta, files -> files },
        ch_bcftools_stats.stats.map { _meta, files -> files },
        gatk_out.mosdepth_summary.map { _meta, files -> files },
        gatk_out.mosdepth_global.map { _meta, files -> files },
        gatk_out.insert_size_metrics.map { _meta, files -> files }
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'pathogenepidemiology_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'pathogenepidemiology'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [ report ] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
