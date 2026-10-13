process GUNZIP_REFFASTA {
    input:
    path ref_gz
    
    output:
    path "${ref_gz.baseName}", emit: gzipped
    
    script:
    """
    gunzip -c ${ref_gz} > ${ref_gz.baseName}
    """
}