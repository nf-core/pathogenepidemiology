include { WGET as DOWNLOAD_QUERYFASTA } from '../../../modules/nf-core/wget/main'
include { WGET as DOWNLOAD_HOSTFASTA  } from '../../../modules/nf-core/wget/main'
include { GUNZIP as GUNZIP_QUERYFASTA } from '../../../modules/nf-core/gunzip/main'
include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'

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
    // TODO: Avoid query/host publication filename collisions.

    // A reference is shared by every sample and interval: keep reusable value channels.
    ch_query = DOWNLOAD_QUERYFASTA(channel.value(query)).outfile
    if (hosturl) {
        hostfasta = DOWNLOAD_HOSTFASTA(channel.value(referenceUriToTuple(hosturl))).outfile.map { meta, fasta -> fasta }
    } else {
        hostfasta = channel.value([])
    }
    queryfasta = GUNZIP_QUERYFASTA(ch_query).gunzip
        .map { meta, fasta -> tuple([id: fasta.baseName], fasta) }
    queryfai = SAMTOOLS_FAIDX(queryfasta.map { meta, fasta -> tuple(meta, fasta, []) }, false).fai
    querydict = GATK4_CREATESEQUENCEDICTIONARY(queryfasta).dict

    emit:
    queryfasta
    queryfai
    querydict
    hostfasta
}
