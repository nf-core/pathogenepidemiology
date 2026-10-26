include { WGET as DOWNLOAD_QUERYFASTA } from '../../modules/nf-core/wget/main'
include { WGET as DOWNLOAD_HOSTFASTA  } from '../../modules/nf-core/wget/main'
include { GUNZIP as GUNZIP_QUERYFASTA } from '../../modules/nf-core/gunzip/main'
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main'

// WGET takes metadata, an explicit URL (part of the task cache key), and a suffix.
def downloadReference(url) {
    def reference = file(new URI(url.toString()).path)
    return tuple([id: reference.baseName], url, reference.extension)
}

workflow PREPARE_REFERENCES {
    take:
    queryurl
    hosturl

    main:
    def query = downloadReference(queryurl)
    def host = downloadReference(hosturl)
    // All references are published together: reject collisions before launching tasks.
    def query_name = "${query[0].id}.${query[2]}"
    def host_name = "${host[0].id}.${host[2]}"
    def query_plain = file(query_name).baseName
    if (host_name in [query_name, query_plain, "${query_plain}.fai"]) {
        error 'Query and host reference filenames conflict in <outdir>/download; use URLs with distinct basenames.'
    }

    ch_query = DOWNLOAD_QUERYFASTA(channel.of(query)).outfile
    hostfasta = DOWNLOAD_HOSTFASTA(channel.of(host)).outfile.map { meta, fasta -> fasta }
    queryfasta = GUNZIP_QUERYFASTA(ch_query).gunzip
        .map { meta, fasta -> tuple([id: fasta.baseName], fasta) }
    queryfai = SAMTOOLS_FAIDX(queryfasta.map { meta, fasta -> tuple(meta, fasta, []) }, false).fai

    emit:
    queryfasta
    queryfai
    hostfasta
}
