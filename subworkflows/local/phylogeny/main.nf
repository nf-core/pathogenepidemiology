nextflow.enable.dsl=2

include { BCFTOOLS_CONSENSUS } from '../../../modules/nf-core/bcftools/consensus/main'
include { MAFFT_ALIGN        } from '../../../modules/nf-core/mafft/align/main'
include { IQTREE             } from '../../../modules/nf-core/iqtree/main'
include { SPLIT_CONSENSUS_CONTIGS } from '../../../modules/local/split_consensus_contigs/main'
include { CONCATENATE_ALIGNED_CONTIGS } from '../../../modules/local/concatenate_aligned_contigs/main'

workflow PHYLOGENY {

    take:
        ch_consensus_input
        ch_reference_fai

    main:

    consensus_s = BCFTOOLS_CONSENSUS(ch_consensus_input).fasta

    ch_consensus_cohort = consensus_s
        .map { meta, fasta -> tuple(meta.id, fasta) }
        .collect(flat: false)
        .map { samples ->
            def sorted_samples = samples.sort { a, b -> a[0] <=> b[0] }
            tuple(
                sorted_samples.collect { sample -> sample[0] },
                sorted_samples.collect { sample -> sample[1] }
            )
        }

    ch_reference_fai_path = ch_reference_fai
        .map { _meta, fai -> fai }
        .first()

    split_s = SPLIT_CONSENSUS_CONTIGS(
        ch_consensus_cohort.combine(ch_reference_fai_path)
    )

    ch_contig_alignments = split_s.contigs
        .flatMap { files ->
            files.sort { a, b -> a.name <=> b.name }.collect { fasta ->
                tuple([id: fasta.baseName], fasta)
            }
        }

    mafft_s = MAFFT_ALIGN(
        ch_contig_alignments,
        [[:], []],
        [[:], []],
        [[:], []],
        [[:], []],
        [[:], []],
        false
    ).fas

        ch_sample_ids = ch_consensus_cohort
        .map { sample_ids, _fastas -> [sample_ids: sample_ids] }

    ch_aligned_contigs = mafft_s
        .map { _meta, alignment -> alignment }
        .collect()
        .map { alignments -> [alignments: alignments] }

    concatenate_input = ch_sample_ids
        .combine(ch_aligned_contigs)
        .combine(ch_reference_fai_path)
        .map { sample_data, alignment_data, reference_fai ->
            [sample_data.sample_ids, alignment_data.alignments, reference_fai]
        }

    concatenate_s = CONCATENATE_ALIGNED_CONTIGS(concatenate_input)

    iqtree_s = IQTREE(
        concatenate_s.alignment.map { alignment -> tuple([id: 'cohort'], alignment, []) },
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        []
    ).phylogeny

    emit:
    phylogeny = iqtree_s
}
