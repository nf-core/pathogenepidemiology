include { WGET as DOWNLOAD_QUERYFASTA } from '../../../modules/nf-core/wget/main'
include { WGET as DOWNLOAD_HOSTFASTA  } from '../../../modules/nf-core/wget/main'
include { GUNZIP as GUNZIP_QUERYFASTA } from '../../../modules/nf-core/gunzip/main'
include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx/main'

// WGET takes metadata, an explicit URL (part of the task cache key), and a suffix.
def referenceUriToTuple(url) {
    def path = file(new URI(url).path)
    tuple([id: path.baseName], url, path.extension)
}

workflow PREPARE_REFERENCES {
    take:
    queryurl
    hosturl

    main:
    def query = referenceUriToTuple(queryurl)
    def host = referenceUriToTuple(hosturl)
    // TODO: Avoid query/host publication filename collisions.

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
