nextflow.enable.dsl=2

include { BCFTOOLS_CONSENSUS } from '../../../modules/nf-core/bcftools/consensus/main'
include { MAFFT_ALIGN        } from '../../../modules/nf-core/mafft/align/main'
include { IQTREE             } from '../../../modules/nf-core/iqtree/main'

workflow PHYLOGENY {

    take:
    ch_consensus_input

    main:

    consensus_s = BCFTOOLS_CONSENSUS(ch_consensus_input).fasta

    ch_consensus_multi = consensus_s
        .map { meta, fasta -> fasta }
        .collectFile(name: 'consensus.fasta')

    mafft_s = MAFFT_ALIGN(
        ch_consensus_multi.map { fasta -> tuple([id: 'cohort'], fasta) },
        [[:], []],
        [[:], []],
        [[:], []],
        [[:], []],
        [[:], []],
        false
    ).fas

    iqtree_s = IQTREE(
        mafft_s.map { meta, alignment -> tuple(meta, alignment, []) },
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
