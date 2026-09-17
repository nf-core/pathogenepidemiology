process GATK4_MARKDUPLICATES_LOCAL {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e3/e3d753d93f57969fe76b8628a8dfcd23ef44bccd08c4ced7089c1f94bf47c89f/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel_htslib_samtools:d3becb6465454c35'}"

    input:
    tuple val(meta), path(bam)
    path fasta
    path fasta_fai
    path core_bed              // NEW

    output:
    tuple val(meta), path("*.pf.bam"),     emit: bam
    tuple val(meta), path("*.pf.bam.bai"), emit: bai
    tuple val(meta), path("*.metrics"),    emit: metrics
    tuple val("${task.process}"), val('gatk4'),    eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"),      topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args  ?: ''
    def args2  = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.sorted.dup"
    def prefix_bam = "${prefix}.bam"
    def prefix_pf  = "${prefix}.pf.bam"

    def input_list = bam.collect { bam_ -> "--INPUT ${bam_}" }.join(' ')
    def reference  = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        MarkDuplicates \\
        ${input_list} \\
        --OUTPUT ${prefix_bam} \\
        --METRICS_FILE ${prefix}.metrics \\
        --TMP_DIR . \\
        ${reference} \\
        ${args}

    samtools view -b ${args2} -L ${core_bed} -o ${prefix_pf} ${prefix_bam}
    samtools index ${prefix_pf}

    rm ${prefix_bam}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.sorted.dup"
    """
    touch ${prefix}.pf.bam
    touch ${prefix}.pf.bam.bai
    touch ${prefix}.metrics
    """
}