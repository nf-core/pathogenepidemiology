process BBDUK_CUSTOM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data' :
        'community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9' }"

    input:
    tuple val(meta), path(reads)
    path references

    output:
    tuple val(meta), path("*_${task.ext.suffix ?: 'trimmed'}.fastq.gz"), emit: reads
    tuple val(meta), path("*.${task.ext.suffix ?: 'trimmed'}.bbduk.log"), emit: log
    tuple val(meta), path("*.${task.ext.suffix ?: 'trimmed'}.discarded.fastq.gz"), emit: discarded
    tuple val(meta), path("*.${task.ext.suffix ?: 'trimmed'}.stats.txt"), emit: stats
    tuple val("${task.process}"), val('bbmap'), eval('bbversion.sh | grep -v "Duplicate cpuset"'), emit: versions_bbmap, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: 'trimmed'
    def raw      = meta.single_end ? "in=${reads[0]}" : "in1=${reads[0]} in2=${reads[1]}"
    def output_reads = meta.single_end ? "out=${prefix}_${suffix}.fastq.gz" : "out1=${prefix}_1_${suffix}.fastq.gz out2=${prefix}_2_${suffix}.fastq.gz"
    def refs = [task.ext.reference ?: [], references].flatten().findAll { ref -> ref }.join(',')
    def reference_arg = refs ? "ref=${refs}" : ''
    // outm holds trimming/length rejects in trimming mode, whole-read matches in filtering mode.
    // TODO: Evaluate upstream nf-core bbmap/bbduk support for stats and discarded outputs.
    def discarded = meta.single_end ? "outm=${prefix}.${suffix}.discarded.fastq.gz" : "outm1=${prefix}_1.${suffix}.discarded.fastq.gz outm2=${prefix}_2.${suffix}.discarded.fastq.gz"
    """
    bbduk.sh \\
        -Xmx${task.memory.toGiga()}g \\
        tossbrokenreads=t qin=33 \\
        $raw \\
        $output_reads \\
        $discarded \\
        stats=${prefix}.${suffix}.stats.txt \\
        $reference_arg \\
        $args \\
        threads=$task.cpus \\
        &> ${prefix}.${suffix}.bbduk.log
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: 'trimmed'
    def output_command = meta.single_end ?
        "echo '' | gzip > ${prefix}_${suffix}.fastq.gz ; echo '' | gzip > ${prefix}.${suffix}.discarded.fastq.gz" :
        "echo '' | gzip > ${prefix}_1_${suffix}.fastq.gz ; echo '' | gzip > ${prefix}_2_${suffix}.fastq.gz ; echo '' | gzip > ${prefix}_1.${suffix}.discarded.fastq.gz ; echo '' | gzip > ${prefix}_2.${suffix}.discarded.fastq.gz"
    """
    touch ${prefix}.${suffix}.bbduk.log
    touch ${prefix}.${suffix}.stats.txt
    $output_command
    """
}
