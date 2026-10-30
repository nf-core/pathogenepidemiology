// Copy selected inputs into task outputs: publication is never a downstream input.
process REUSE_REFERENCE {
    tag "${meta.id}"
    label 'process_single'

    input:
    // Keep the original basename separate from Nextflow's staged input/ path.
    tuple val(meta), path(asset, stageAs: 'input/*'), val(filename)

    output:
    tuple val(meta), path(filename), emit: reference

    script:
    """
    cp ${asset} ${filename}
    """
}

process REUSE_BWAMEM3 {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(assets, stageAs: 'input/*')

    output:
    tuple val(meta), path('bwamem3'), emit: index

    script:
    """
    mkdir bwamem3
    cp input/* bwamem3/
    """
}
