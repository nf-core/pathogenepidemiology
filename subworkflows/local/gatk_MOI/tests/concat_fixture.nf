// Synthetic records only: one site per chromosome and numeric header order.
process CONCAT_FIXTURE {
    container 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'
    cpus 1
    memory 1.GB

    input:
    tuple val(meta), val(ready)

    output:
    tuple val(meta), path('chr*.vcf.gz'), emit: vcfs
    tuple val(meta), path('chr*.vcf.gz.tbi'), emit: indexes

    script:
    """
    for n in \$(seq 1 14); do
        {
            printf '##fileformat=VCFv4.2\\n'
            for c in \$(seq 1 14); do
                printf '##contig=<ID=chr%s,length=1000>\\n' "\$c"
            done
            printf '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n'
            printf 'chr%s\\t10\\t.\\tA\\tC\\t30\\tPASS\\t.\\n' "\$n"
        } | bcftools view -Oz --write-index=tbi -o "chr\$n.vcf.gz"
    done
    """
}
