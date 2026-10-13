process DOWNLOAD_REFFASTA {
    storeDir "${params.outdir}/download"
    
    input:
    val ref_url
    
    output:
    path "${ref_url.toString().split('/').last()}", emit:fasta
    
    script:
    """
    wget -O ${ref_url.toString().split('/').last()} ${ref_url}
    """
}