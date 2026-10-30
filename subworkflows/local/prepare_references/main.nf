include { WGET as DOWNLOAD_QUERYFASTA } from '../../../modules/nf-core/wget/main'
include { WGET as DOWNLOAD_HOSTFASTA  } from '../../../modules/nf-core/wget/main'
include { GUNZIP as GUNZIP_QUERYFASTA } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_HOSTFASTA  } from '../../../modules/nf-core/gunzip/main'
include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { MINIMAP2_INDEX } from '../../../modules/nf-core/minimap2/index/main'
include { BWAMEM3_INDEX } from '../../../modules/nf-core/bwamem3/index/main'

def adjacentReference(reference, name) {
    def asset = reference.resolveSibling(name)
    // A present but invalid companion is an input error, not a cache miss.
    if (java.nio.file.Files.exists(asset, java.nio.file.LinkOption.NOFOLLOW_LINKS)) {
        if (!java.nio.file.Files.isRegularFile(asset) || !java.nio.file.Files.isReadable(asset)) {
            error "Reference companion must be a regular readable file: ${asset}"
        }
        return asset
    }
    return null
}

def localReference(value, option) {
    if (!(value instanceof CharSequence) || value.contains('://') || !(value ==~ /.+\.(fa|fasta|fna)(\.gz)?/)) {
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
    def companions = [:]
    if (local_query != null) {
        def reference = local_query[1]
        def plain = reference.name.replaceFirst(/\.gz$/, '')
        def stem = plain.replaceFirst(/\.(fa|fasta|fna)$/, '')
        companions.fai = adjacentReference(reference, "${plain}.fai")
        companions.dict = adjacentReference(reference, "${stem}.dict")
        companions.mmi = adjacentReference(reference, "${stem}.mmi")
        def bwa = ['0123', 'amb', 'ann', 'bwt.2bit.64', 'pac'].collect {
            adjacentReference(reference, "${plain}.${it}")
        }
        // BWA indexes are atomic: never mix supplied and rebuilt sidecars.
        companions.bwa = bwa.every { it != null } ? bwa : null
    }
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
    queryfai = companions.fai != null
        ? queryfasta.map { meta, fasta -> tuple(meta, companions.fai) }
        : SAMTOOLS_FAIDX(queryfasta.map { meta, fasta -> tuple(meta, fasta, []) }, false).fai
    querydict = companions.dict != null
        ? queryfasta.map { meta, fasta -> tuple(meta, companions.dict) }
        : GATK4_CREATESEQUENCEDICTIONARY(queryfasta).dict
    querymmi = companions.mmi != null
        ? queryfasta.map { meta, fasta -> tuple(meta, companions.mmi) }
        : MINIMAP2_INDEX(queryfasta).index
    // BWAMEM3_MEM accepts path collections as well as generated index directories.
    // Stage only the exact five companions, never their parent with unrelated indexes.
    querybwamem3 = companions.bwa != null
        ? queryfasta.map { meta, fasta -> tuple(meta, companions.bwa) }
        : BWAMEM3_INDEX(queryfasta).index

    emit:
    queryfasta
    queryfai
    querydict
    querymmi
    querybwamem3
    hostfasta
}
