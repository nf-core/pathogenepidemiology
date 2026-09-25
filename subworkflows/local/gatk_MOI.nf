// GATK4 short-read variant-calling subworkflow - refactor from Niare et al
//
// STATUS: draft scaffold -- wires the full chain (CleanSam -> ... -> ApplyVQSR)


//include { SAMTOOLS_INDEX                  } from '../../modules/nf-core/samtools/index/main'
include { GATK4_CLEANSAM                  } from '../../modules/nf-core/gatk4/cleansam/main'
include { PICARD_SORTSAM                  } from '../../modules/nf-core/picard/sortsam/main'
include { GATK4_CREATESEQUENCEDICTIONARY  } from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_MARKDUPLICATES_CUSTOM      } from '../../modules/local/gatk4/markduplicates_custom/main'
include { MOSDEPTH                        } from '../../modules/nf-core/mosdepth/main'
include { PICARD_COLLECTINSERTSIZEMETRICS } from '../../modules/nf-core/picard/collectinsertsizemetrics/main'
include { GATK4_HAPLOTYPECALLER           } from '../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_GENOMICSDBIMPORT          } from '../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS             } from '../../modules/nf-core/gatk4/genotypegvcfs/main'
// step: gatk GatherVCFs has no module. bcftools concat serves same purpose
include { BCFTOOLS_CONCAT } from '../../modules/nf-core/bcftools/concat/main'
include { 
  GATK4_VARIANTRECALIBRATOR as GATK4_VARRECAL_INDELS      
          } from '../../modules/nf-core/gatk4/variantrecalibrator/main'
include { 
  GATK4_APPLYVQSR as GATK4_VQSR_INDELS
  } from '../../modules/nf-core/gatk4/applyvqsr/main'
include { 
  GATK4_VARIANTRECALIBRATOR as GATK4_VARRECAL_SNPS      
          } from '../../modules/nf-core/gatk4/variantrecalibrator/main'
include { 
  GATK4_APPLYVQSR as GATK4_VQSR_SNPS
  } from '../../modules/nf-core/gatk4/applyvqsr/main'



workflow GATK_MOI {

    take:
    ch_aligned_s          // channel: [ meta, bam ]         sorted short-read BAM from BWAMEM3_MEM (sort=true)
    ch_queryfasta         // channel: [ meta, fasta ]       from PREPARE_REFERENCES.out.queryfasta
    ch_queryfai           // channel: [ meta, fai ]         from PREPARE_REFERENCES.out.queryfai
    ch_vqsr_resource      // channel: path(training_vcf)    equivalent of the paper's Strains.vcf.gz
    ch_vqsr_resource_tbi  // channel: path(training_vcf_tbi)


    main:

    // Prepare channels for specific use cases:

    ch_versions = Channel.empty() // accumulates tool versions, which can be emitted as 
                                  // "versions_gatk4", "versions_gatk" or "versions_picard"

    GATK4_CREATESEQUENCEDICTIONARY(ch_queryfasta)
    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions_gatk4)

    ch_querydict = GATK4_CREATESEQUENCEDICTIONARY.out.dict

    // Create val versions (metadata-free) of channels for modules that expect bare paths
    ch_fasta_val = ch_queryfasta.map { meta, fasta -> fasta }.first()
    ch_fai_val     = ch_queryfai.map { meta, fai -> fai }.first()
    ch_dict_val   = ch_querydict.map { meta, dict -> dict }.first()

    // channel for Pf chromosomal intervals
    // TODO: launchDir -> baseDir when this is exec'd through pipeline main script
    ch_intervals = Channel.fromPath("${launchDir}/assets/intervals/core_chr*.list") 
        .map { list ->
            def chr_id = list.baseName.replaceAll(/core_chr0*/, 'chr') // e.g., "core_chr01" -> "chr1"
            tuple([id: chr_id], list)
        }

    // Start pipeline proper:

    // QC steps

    GATK4_CLEANSAM(
        ch_aligned_s,
        ch_queryfasta.join(ch_queryfai, by: 0).first()   // tuple(meta2, fasta, fasta_index) -- module needs both together
    )
    ch_versions = ch_versions.mix(GATK4_CLEANSAM.out.versions_gatk)

    PICARD_SORTSAM(
        GATK4_CLEANSAM.out.bam,
        'coordinate'
    )
    ch_versions = ch_versions.mix(PICARD_SORTSAM.out.versions_picard)

    // TODO: change to projectDir or baseDir before pipeline is shipped
    //       /when it should run from pipelineroot/main.nf
    def core_bed = file("${launchDir}/assets/Pf3D7_core.bed")

    GATK4_MARKDUPLICATES_CUSTOM(
    PICARD_SORTSAM.out.bam,
    ch_fasta_val,
    ch_fai_val,
    core_bed
    )
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES_CUSTOM.out.versions_gatk4)
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES_CUSTOM.out.versions_samtools)

    ch_dedup_bam = GATK4_MARKDUPLICATES_CUSTOM.out.bam
        .join(GATK4_MARKDUPLICATES_CUSTOM.out.bai, by: 0)


    PICARD_COLLECTINSERTSIZEMETRICS(
      ch_dedup_bam.map { meta, bam, bai -> tuple(meta, bam) }
      )
    ch_versions = ch_versions.mix(PICARD_COLLECTINSERTSIZEMETRICS.out.versions_picard)

    // NOTE: paper's DepthOfCoverage has no nf-core module. Using mosdepth instead.

    ch_mosdepth_input = ch_dedup_bam
    .map { meta, bam, bai -> tuple(meta, bam, bai, core_bed) }

    MOSDEPTH(
        ch_mosdepth_input,   // [meta, bam, bai, bed]
        ch_queryfasta.first(),            // [meta2, fasta]
        []                   // quantize_labels — empty means no quantization
    )
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions_mosdepth)
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions_gzip)

    ch_hc_input = ch_dedup_bam
        .combine(ch_intervals)
        .map { meta, bam, bai, interval_meta, interval_file ->
            def combined_meta = meta + [chr: interval_meta.id]
            tuple(combined_meta, bam, bai, interval_file, [])
        }


    GATK4_HAPLOTYPECALLER(
    ch_hc_input,
    ch_queryfasta.first(),
    ch_queryfai.first(),
    ch_querydict.first(),
    [ [:], [] ],
    [ [:], [] ]
    )
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions_gatk4)


    // Regroup gVCFs so each item holds all samples for one chromosome
    ch_gdbi_input = GATK4_HAPLOTYPECALLER.out.vcf
    .join(GATK4_HAPLOTYPECALLER.out.tbi, by: 0)
    .map { meta, vcf, tbi -> tuple([chr: meta.chr], vcf, tbi) }
    .groupTuple(by: 0)
    .map { chr_meta, vcfs, tbis ->
        def chr_num  = chr_meta.chr.replace('chr', '')
        def interval = file("${launchDir}/assets/intervals/core_chr${chr_num}.list")
        tuple(
            [id: chr_meta.chr],   // ← meta needs an id, not just chr
            vcfs,                 // ← list, which module iterates over. produces one --variant flag per file.
            tbis,
            interval,             // ← real .list file, not []
            "",                   // ← interval_value slot, empty
            []                    // ← wspace slot, empty
        )
    }

    GATK4_GENOMICSDBIMPORT(
        ch_gdbi_input,
        false,  // run_intlist
        false,  // run_updatewspace
        false    // input_map
    )
    ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions_gatk4)

    GATK4_GENOTYPEGVCFS(
        GATK4_GENOMICSDBIMPORT.out.genomicsdb.map { meta, db -> tuple(meta, db, [], [], []) },
        ch_queryfasta.first(),
        ch_queryfai.first(),
        ch_querydict.first(),
        [ [:], [] ],
        [ [:], [] ]
    )
    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions_gatk4)

    

    ch_raw_vcf = GATK4_GENOTYPEGVCFS.out.vcf.join(GATK4_GENOTYPEGVCFS.out.tbi, by: 0)


    GATK4_VARRECAL_INDELS(
        ch_raw_vcf,
        ch_vqsr_resource,
        ch_vqsr_resource_tbi,
        [ '-resource:Strains,known=true,training=true,truth=true,prior=15.0 Strains.2kb.vcf.gz' ],
        ch_fasta_val,
        ch_fai_val,
        ch_dict_val
    )
    ch_versions = ch_versions.mix(GATK4_VARRECAL_INDELS.out.versions_gatk4)

    GATK4_VQSR_INDELS(
        ch_raw_vcf
            .join(GATK4_VARRECAL_INDELS.out.recal, by: 0)
            .join(GATK4_VARRECAL_INDELS.out.idx, by: 0)
            .join(GATK4_VARRECAL_INDELS.out.tranches, by: 0),
        ch_fasta_val,
        ch_fai_val,
        ch_dict_val
    )
    ch_versions = ch_versions.mix(GATK4_VQSR_INDELS.out.versions_gatk4)

    ch_indelrecal_vcf = GATK4_VQSR_INDELS.out.vcf.join(GATK4_VQSR_INDELS.out.tbi, by: 0)

    GATK4_VARRECAL_SNPS(
        ch_indelrecal_vcf,
        ch_vqsr_resource,
        ch_vqsr_resource_tbi,
        [ '-resource:Strains,known=true,training=true,truth=true,prior=15.0 Strains.2kb.vcf.gz' ],
        ch_fasta_val,
        ch_fai_val,
        ch_dict_val
    )
    ch_versions = ch_versions.mix(GATK4_VARRECAL_SNPS.out.versions_gatk4)

    GATK4_VQSR_SNPS(
        ch_indelrecal_vcf
            .join(GATK4_VARRECAL_SNPS.out.recal, by: 0)
            .join(GATK4_VARRECAL_SNPS.out.idx, by: 0)
            .join(GATK4_VARRECAL_SNPS.out.tranches, by: 0),
        ch_fasta_val,
        ch_fai_val,
        ch_dict_val
    )
    ch_versions = ch_versions.mix(GATK4_VQSR_SNPS.out.versions_gatk4)

    ch_vqsr_output = GATK4_VQSR_SNPS.out.vcf.join(GATK4_VQSR_SNPS.out.tbi, by: 0)


    // Concatenate the 14 per-chromosome cohort VCFs into one
    ch_concat_input = ch_vqsr_output
        .map { meta, vcf, tbi -> tuple(meta.id, vcf, tbi) }
        .toSortedList { a, b -> (a[0] - 'chr') as Integer <=> (b[0] - 'chr') as Integer }
        .map { items ->
            def vcfs = items.collect { it[1] }
            def tbis = items.collect { it[2] }
            tuple([id: 'cohort'], vcfs, tbis)
        }

    BCFTOOLS_CONCAT(ch_concat_input)
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions_bcftools)

    varcalls_s = BCFTOOLS_CONCAT.out.vcf.join(BCFTOOLS_CONCAT.out.index, by: 0)

    emit:
    varcalls_s          = varcalls_s
    insert_size_metrics = PICARD_COLLECTINSERTSIZEMETRICS.out.metrics
    versions            = ch_versions
    mosdepth_summary    = MOSDEPTH.out.summary_txt
    mosdepth_global     = MOSDEPTH.out.global_txt
}
/*
// for standalone testing while we work on Issue #13
include { PREPARE_REFERENCES } from './prepare_references'

workflow {
    ch_aligned_s = Channel.fromPath("${launchDir}/results/aligned_s/*.bam")
        .map { bam -> tuple([id: bam.baseName], bam) }

    PREPARE_REFERENCES(params.queryurl, params.hosturl)

    gatk_out = GATK_MOI(
        ch_aligned_s,
        PREPARE_REFERENCES.out.queryfasta,
        PREPARE_REFERENCES.out.queryfai,
        Channel.fromPath("${launchDir}/assets/Strains.2kb.vcf.gz"),   
        Channel.fromPath("${launchDir}/assets/Strains.2kb.vcf.gz.tbi")  
    )

    gatk_out.varcalls_s.view { meta, vcf -> "gVCF: ${meta.id} -> ${vcf}" }
}
*/