include { WGET as DOWNLOAD_QUERYFASTA } from '../../../modules/nf-core/wget/main'
include { WGET as DOWNLOAD_HOSTFASTA  } from '../../../modules/nf-core/wget/main'
include { GUNZIP as GUNZIP_QUERYFASTA } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_HOSTFASTA  } from '../../../modules/nf-core/gunzip/main'
include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'

def localReference(value, option) {
    if (!(value instanceof CharSequence) || !(value ==~ /.+\.(fa|fasta|fna)(\.gz)?/)) {
        error "--${option} must be a single local FASTA file (.fa, .fasta or .fna, optionally .gz)."
    }
    def reference = file(value, checkIfExists: true, glob: false)
    if (!reference.isFile()) {
        error "--${option} must be a single file: ${value}"
    }
    // Vendored processes interpolate staged basenames into shell commands.
    if (!(reference.name ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) {
        error "--${option} basename must start with a letter, digit or underscore and contain only letters, digits, dots, underscores and hyphens: ${reference.name}"
    }
    return tuple([id: reference.baseName], reference)
}

// WGET takes metadata, an explicit URL (part of the task cache key), and a suffix.
def referenceUriToTuple(url) {
    def path = file(file(new URI(url.toString()).path).name)
    tuple([id: path.baseName], url, path.extension)
}

workflow PREPARE_REFERENCES {
    take:
    queryurl
    hosturl
    query_fasta
    host_fasta

    main:
    // Explicit local files take precedence independently; validate before launching tasks.
    def local_query = query_fasta != null ? localReference(query_fasta, 'query_fasta') : null
    def local_host = host_fasta != null ? localReference(host_fasta, 'host_fasta') : null
    def query = local_query != null ? local_query : referenceUriToTuple(queryurl)
    def host = local_host != null ? local_host : (hosturl ? referenceUriToTuple(hosturl) : null)
    // TODO: Avoid query/host publication filename collisions.

    // A reference is shared by every sample and interval: keep reusable value channels.
    ch_query = local_query != null ? channel.value(query) : DOWNLOAD_QUERYFASTA(channel.value(query)).outfile
    if (host != null) {
        ch_host = local_host != null ? channel.value(host) : DOWNLOAD_HOSTFASTA(channel.value(host)).outfile
        def host_gzip = local_host != null ? host[1].extension == 'gz' : host[2] == 'gz'
        hostfasta = (host_gzip ? GUNZIP_HOSTFASTA(ch_host).gunzip : ch_host)
            .map { meta, fasta -> fasta }
    } else {
        hostfasta = channel.value([])
    }
    def query_gzip = local_query != null ? query[1].extension == 'gz' : query[2] == 'gz'
    queryfasta = (query_gzip ? GUNZIP_QUERYFASTA(ch_query).gunzip : ch_query)
        .map { meta, fasta -> tuple([id: fasta.baseName], fasta) }
    queryfai = SAMTOOLS_FAIDX(queryfasta.map { meta, fasta -> tuple(meta, fasta, []) }, false).fai
    querydict = GATK4_CREATESEQUENCEDICTIONARY(queryfasta).dict

    emit:
    queryfasta
    queryfai
    querydict
    hostfasta
}
