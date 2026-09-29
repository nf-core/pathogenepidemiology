// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process HMMIBDRS {
    tag "$meta.id"
    label 'process_high'

    //conda "${moduleDir}/environment.yml" // TODO: reintroduce conda directive when the package is accepted
    // TODO: rewrite this statement in nf-core house style when there is a Sequera-hosted alternative
    container "docker://bguo068/hmmibd-rs:v0.1.5"

    input:
    tuple val(meta), path(genotypes), path(genotypes_index)

    output:
    tuple val(meta), path("*.hmm.txt"),       emit: hmm
    tuple val(meta), path("*.hmm_final.txt"), emit: hmm_final
    tuple val("${task.process}"), val('hmmibd-rs'), eval("hmmibd-rs --version | sed 's/^hmmibd-rs //'"), topic: versions, emit: versions_hmmibdrs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // For this tool, "--num-threads 0" is the default and uses all threads
    """
    hmmibd-rs \\
        --data-file1 ${genotypes} \\
        --output ${prefix} \\
        --num-threads ${task.cpus} \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/624977dfaf562211e68a8a868ca80acc8461f1ac/modules/nf-core/cutadapt/main.nf#L34-L46
    //               Complex example: https://github.com/nf-core/modules/blob/88d43dad73a675e66bff49ebb57fe657a5909018/modules/nf-core/bedtools/split/main.nf#L32-L43
    """
    echo $args

    touch ${prefix}.hmm.txt
    touch ${prefix}.hmm_final.txt
    """
}