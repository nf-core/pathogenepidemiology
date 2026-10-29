include { CONCAT_COHORT } from '../main.nf'

process MAKE_VCFS {
    label 'process_single'
    container 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'

    output:
    tuple path('chr2.vcf.gz'), path('chr2.vcf.gz.tbi'), emit: chr2
    tuple path('chr10.indel.snp.recal.vcf.gz'), path('chr10.indel.snp.recal.vcf.gz.tbi'), emit: chr10
    tuple path('chr1.vcf.gz'), path('chr1.vcf.gz.tbi'), emit: chr1
    path 'query.fa.fai', emit: fai

    script:
    """
    printf 'chr2\\t1000\\t0\\t60\\t61\\nchr10\\t1000\\t0\\t60\\t61\\nchr1\\t1000\\t0\\t60\\t61\\n' > query.fa.fai
    cat > header <<'HEADER'
    ##fileformat=VCFv4.2
    ##contig=<ID=chr2,length=1000>
    ##contig=<ID=chr10,length=1000>
    ##contig=<ID=chr1,length=1000>
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2
    HEADER
    for contig in chr2 chr10 chr1; do
        name=\$contig
        if [[ \$contig == chr10 ]]; then name=chr10.indel.snp.recal; fi
        { cat header; printf '%s\\t10\\t.\\tA\\tC\\t60\\tPASS\\t.\\tGT\\t0/1\\t1/1\\n' "\$contig"; } | bcftools view -Oz -o \$name.vcf.gz
        bcftools index -t \$name.vcf.gz
    done
    """
}

process CHECK_COHORT {
    label 'process_single'
    container 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'

    input:
    tuple val(meta), path(vcf), path(index)

    output:
    path 'records.txt', emit: records
    path 'samples.txt', emit: samples
    path 'indexed.txt', emit: indexed

    script:
    """
    bcftools query -f '%CHROM\t%POS\t[%GT\t]\n' ${vcf} > records.txt
    bcftools query -l ${vcf} > samples.txt
    bcftools index -n '${vcf}##idx##${index}' > indexed.txt
    """
}

workflow TEST_CONCAT_COHORT {
    main:
    MAKE_VCFS()
    // Deliberately mix completion order, lexical order and VQSR/bypass suffixes.
    ch_chr10 = MAKE_VCFS.out.chr10.map { vcf, index -> tuple([id: 'chr10'], vcf, index) }
    ch_chr1 = MAKE_VCFS.out.chr1.map { vcf, index -> tuple([id: 'chr1'], vcf, index) }
    ch_chr2 = MAKE_VCFS.out.chr2.map { vcf, index -> tuple([id: 'chr2'], vcf, index) }
    CONCAT_COHORT(ch_chr10.mix(ch_chr1, ch_chr2), MAKE_VCFS.out.fai.map { fai -> tuple([id: 'query'], fai) })
    CHECK_COHORT(CONCAT_COHORT.out.vcf.join(CONCAT_COHORT.out.index, by: 0))

    emit:
    records = CHECK_COHORT.out.records
    samples = CHECK_COHORT.out.samples
    indexed = CHECK_COHORT.out.indexed
}
