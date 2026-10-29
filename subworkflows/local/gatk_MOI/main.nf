// GATK4 short-read variant-calling subworkflow - refactor from Niare et al
//
// STATUS: draft scaffold -- wires the full chain (CleanSam -> ... -> ApplyVQSR)


include { GATK4_CLEANSAM                                     } from '../../../modules/nf-core/gatk4/cleansam/main'
include { PICARD_SORTSAM                                     } from '../../../modules/nf-core/picard/sortsam/main'
include { GATK4_CREATESEQUENCEDICTIONARY                     } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_MARKDUPLICATES_CUSTOM                        } from '../../../modules/local/gatk4/markduplicates_custom/main'
include { MOSDEPTH                                           } from '../../../modules/nf-core/mosdepth/main'
include { PICARD_COLLECTINSERTSIZEMETRICS                    } from '../../../modules/nf-core/picard/collectinsertsizemetrics/main'
include { GATK4_HAPLOTYPECALLER                              } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_GENOMICSDBIMPORT                             } from '../../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS                                } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'
include { BCFTOOLS_CONCAT                                    } from '../../../modules/nf-core/bcftools/concat/main'
include {
    GATK4_VARIANTRECALIBRATOR as GATK4_VARRECAL_INDELS
} from '../../../modules/nf-core/gatk4/variantrecalibrator/main'
include {
    GATK4_APPLYVQSR as GATK4_VQSR_INDELS
} from '../../../modules/nf-core/gatk4/applyvqsr/main'
include {
    GATK4_VARIANTRECALIBRATOR as GATK4_VARRECAL_SNPS
} from '../../../modules/nf-core/gatk4/variantrecalibrator/main'
include {
    GATK4_APPLYVQSR as GATK4_VQSR_SNPS
} from '../../../modules/nf-core/gatk4/applyvqsr/main'


workflow GATK_MOI {
    take:
    ch_aligned_s
    ch_queryfasta
    ch_queryfai
    ch_vqsr_resource
    ch_vqsr_resource_tbi
    ch_core_bed
    ch_intervals

    main:

    GATK4_CREATESEQUENCEDICTIONARY(ch_queryfasta)

    ch_querydict = GATK4_CREATESEQUENCEDICTIONARY.out.dict

    ch_fasta_val = ch_queryfasta.map { _meta, fasta -> fasta }
    ch_fai_val = ch_queryfai.map { _meta, fai -> fai }
    ch_dict_val = ch_querydict.map { _meta, dict -> dict }

    // Both references are values; map preserves reuse without a queue-producing join.
    ch_reference = ch_queryfasta.map { meta, fasta -> tuple(meta, fasta, ch_queryfai.val[1]) }



    GATK4_CLEANSAM(
        ch_aligned_s,
        ch_reference,
    )

    PICARD_SORTSAM(GATK4_CLEANSAM.out.bam, 'coordinate')

    GATK4_MARKDUPLICATES_CUSTOM(
        PICARD_SORTSAM.out.bam,
        ch_fasta_val,
        ch_fai_val,
        ch_core_bed,
    )

    ch_dedup_bam = GATK4_MARKDUPLICATES_CUSTOM.out.bam.join(GATK4_MARKDUPLICATES_CUSTOM.out.bai, by: 0)

    PICARD_COLLECTINSERTSIZEMETRICS(
        ch_dedup_bam.map { meta, bam, _bai -> tuple(meta, bam) }
    )

    ch_mosdepth_input = ch_dedup_bam.combine(ch_core_bed)

    MOSDEPTH(
        ch_mosdepth_input,
        ch_queryfasta,
        [],
    )

    ch_hc_input = ch_dedup_bam
        .combine(ch_intervals)
        .map { meta, bam, bai, interval_meta, interval_file ->
            def combined_meta = meta + [chr: interval_meta.id]
            tuple(combined_meta, bam, bai, interval_file, [])
        }

    GATK4_HAPLOTYPECALLER(
        ch_hc_input,
        ch_queryfasta,
        ch_queryfai,
        ch_querydict,
        [[:], []],
        [[:], []],
    )

    ch_gdbi_input = GATK4_HAPLOTYPECALLER.out.vcf
        .join(GATK4_HAPLOTYPECALLER.out.tbi, by: 0)
        .map { meta, vcf, tbi -> tuple([id: meta.chr], vcf, tbi) }
        .groupTuple(by: 0)
        .join(ch_intervals, by: 0)
        .map { interval_meta, vcfs, tbis, interval ->
            tuple(
                interval_meta,
                vcfs,
                tbis,
                interval,
                "",
                [],
            )
        }

    GATK4_GENOMICSDBIMPORT(
        ch_gdbi_input,
        false,
        false,
        false,
    )

    GATK4_GENOTYPEGVCFS(
        GATK4_GENOMICSDBIMPORT.out.genomicsdb.map { meta, db -> tuple(meta, db, [], [], []) },
        ch_queryfasta,
        ch_queryfai,
        ch_querydict,
        [[:], []],
        [[:], []],
    )

    ch_vcf_with_index = GATK4_GENOTYPEGVCFS.out.vcf
        .join(GATK4_GENOTYPEGVCFS.out.tbi, by: 0)

    // Preserve the existing VQSR size gate, but retain smaller VCFs unrecalibrated.
    ch_raw_vcf = ch_vcf_with_index.filter { _meta, vcf, _tbi -> vcf.size() > 8000 }
    ch_bypass_vcf = ch_vcf_with_index.filter { _meta, vcf, _tbi -> vcf.size() <= 8000 }


    ch_insert_metrics = PICARD_COLLECTINSERTSIZEMETRICS.out.metrics
    ch_mosdepth_sum = MOSDEPTH.out.summary_txt
    ch_mosdepth_glob = MOSDEPTH.out.global_txt


    GATK4_VARRECAL_INDELS(
        ch_raw_vcf,
        ch_vqsr_resource,
        ch_vqsr_resource_tbi,
        ['-resource:Strains,known=true,training=true,truth=true,prior=15.0 Strains.2kb.vcf.gz'],
        ch_fasta_val,
        ch_fai_val,
        ch_dict_val,
    )

    GATK4_VQSR_INDELS(
        ch_raw_vcf.join(GATK4_VARRECAL_INDELS.out.recal, by: 0).join(GATK4_VARRECAL_INDELS.out.idx, by: 0).join(GATK4_VARRECAL_INDELS.out.tranches, by: 0),
        ch_fasta_val,
        ch_fai_val,
        ch_dict_val,
    )

    ch_indelrecal_vcf = GATK4_VQSR_INDELS.out.vcf.join(GATK4_VQSR_INDELS.out.tbi, by: 0)

    GATK4_VARRECAL_SNPS(
        ch_indelrecal_vcf,
        ch_vqsr_resource,
        ch_vqsr_resource_tbi,
        ['-resource:Strains,known=true,training=true,truth=true,prior=15.0 Strains.2kb.vcf.gz'],
        ch_fasta_val,
        ch_fai_val,
        ch_dict_val,
    )

    GATK4_VQSR_SNPS(
        ch_indelrecal_vcf.join(GATK4_VARRECAL_SNPS.out.recal, by: 0).join(GATK4_VARRECAL_SNPS.out.idx, by: 0).join(GATK4_VARRECAL_SNPS.out.tranches, by: 0),
        ch_fasta_val,
        ch_fai_val,
        ch_dict_val,
    )

    ch_vqsr_output = GATK4_VQSR_SNPS.out.vcf.join(GATK4_VQSR_SNPS.out.tbi, by: 0)

    ch_concat_input = ch_vqsr_output
        .mix(ch_bypass_vcf)
        .collect(flat: false)
        .map { items ->
            tuple([id: 'cohort'], items.collect { item -> item[1] }, items.collect { item -> item[2] })
        }

    BCFTOOLS_CONCAT(ch_concat_input)

    varcalls_s = BCFTOOLS_CONCAT.out.vcf.join(BCFTOOLS_CONCAT.out.index, by: 0)

    emit:
    varcalls_s          = varcalls_s
    insert_size_metrics = ch_insert_metrics
    mosdepth_summary    = ch_mosdepth_sum
    mosdepth_global     = ch_mosdepth_glob
}
