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
// step: gatk GatherVCFs
include { 
  GATK4_VARIANTRECALIBRATOR as GATK_VARRECAL_INDELS      
          } from '../../modules/nf-core/gatk4/variantrecalibrator/main'
include { 
  GATK4_APPLYVQSR as GATK4_VQSR_INDELS
  } from '../../modules/nf-core/gatk4/applyvqsr/main'
include { 
  GATK4_VARIANTRECALIBRATOR as GATK_VARRECAL_SNPS      
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


    def core_bed = file("${baseDir}/assets/Pf3D7_core.bed")

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


    // TODO: DECISION: -L arg, either create interval files and channel for chromosomally 
    //       parallelised analysis, or break from original Niare+al pipeline design. 
    GATK4_HAPLOTYPECALLER(
    ch_dedup_bam.map { meta, bam -> tuple(meta, bam, file("${bam}.bai"), [], []) },
    ch_queryfasta.first(),
    ch_queryfai.first(),
    ch_querydict.first(),
    [ [:], [] ],
    [ [:], [] ]
    )
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions_gatk4)

/*
    // ============================================================
    // STAGE D -- joint genotyping across samples
    // paper: GenomicsDBImport (per chromosome) -> GenotypeGVCFs (per
    //        genomic sub-region, run as parallel SLURM jobs) -> GatherVcfs
    // ============================================================

    // TODO real gap: GenomicsDBImport combines gVCFs from ALL samples into
    // one database per interval. The single-sample tuple below is only a
    // placeholder to keep the workflow syntactically chainable -- the real
    // version needs ch_gvcfs grouped/collected across the whole sample set,
    // keyed by chromosome/interval (this is the Nextflow equivalent of the
    // paper's per-chromosome `for i in 1..14` loop -- an interval channel,
    // not a bash loop). Need to design the interval channel first.
    ch_gvcfs_for_import = GATK4_HAPLOTYPECALLER.out.vcf
        .join(GATK4_HAPLOTYPECALLER.out.tbi, by: 0)
        .map { meta, vcf, tbi -> tuple(meta, vcf, tbi, [], [], []) }

    GATK4_GENOMICSDBIMPORT(
        ch_gvcfs_for_import,
        false,  // run_intlist
        false,  // run_updatewspace
        true    // input_map
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

    // NOTE: paper's GatherVcfs (recombine per-region VCFs into one per
    // chromosome) has no nf-core module. BCFTOOLS_CONCAT does the same job
    // and does have one -- worth checking before writing a local module.

    // ============================================================
    // STAGE E -- VQSR filtering: indels first, then SNPs on that output
    // ============================================================

    ch_raw_vcf = GATK4_GENOTYPEGVCFS.out.vcf.join(GATK4_GENOTYPEGVCFS.out.tbi, by: 0)

    // TODO: 'labels' format below is a guess -- confirm against the module's
    // actual expected resource-label syntax before running. Also,
    // ch_vqsr_resource(_tbi) needs a real training VCF -- the paper uses a
    // custom "Strains.vcf.gz"; there isn't an equivalent asset yet
    // so this needs sourcing/creating before Stage E can actually run.
    VARCAL_INDEL(
        ch_raw_vcf,
        ch_vqsr_resource,
        ch_vqsr_resource_tbi,
        [ 'Brown,known=true,training=true,truth=true,prior=15.0' ],
        ch_fasta_val,
        ch_fai_val,
        ch_dict_val
    )

    VQSR_INDEL(
        ch_raw_vcf
            .join(VARCAL_INDEL.out.recal, by: 0)
            .join(VARCAL_INDEL.out.idx, by: 0)
            .join(VARCAL_INDEL.out.tranches, by: 0),
        ch_fasta_val,
        ch_fai_val,
        ch_dict_val
    )

    ch_indelrecal_vcf = VQSR_INDEL.out.vcf.join(VQSR_INDEL.out.tbi, by: 0)

    VARCAL_SNP(
        ch_indelrecal_vcf,
        ch_vqsr_resource,
        ch_vqsr_resource_tbi,
        [ 'Brown,known=true,training=true,truth=true,prior=15.0' ],
        ch_fasta_val,
        ch_fai_val,
        ch_dict_val
    )

    VQSR_SNP(
        ch_indelrecal_vcf
            .join(VARCAL_SNP.out.recal, by: 0)
            .join(VARCAL_SNP.out.idx, by: 0)
            .join(VARCAL_SNP.out.tranches, by: 0),
        ch_fasta_val,
        ch_fai_val,
        ch_dict_val
    )

    varcalls_s = VQSR_SNP.out.vcf.join(VQSR_SNP.out.tbi, by: 0)
*/
    emit:
    //varcalls_s                                                    // channel: [ meta, vcf, tbi ] -- final recalibrated short-read VCF
    varcalls_s = GATK4_HAPLOTYPECALLER.out.vcf
    insert_size_metrics = PICARD_COLLECTINSERTSIZEMETRICS.out.metrics
    versions             = ch_versions
    mosdepth_summary = MOSDEPTH.out.summary_txt
    mosdepth_global  = MOSDEPTH.out.global_txt
}
