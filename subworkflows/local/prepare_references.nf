
include { DOWNLOAD_REFFASTA as DOWNLOAD_QUERYFASTA } from '../../modules/local/download_reffasta/main'
include { DOWNLOAD_REFFASTA as DOWNLOAD_HOSTFASTA  } from '../../modules/local/download_reffasta/main'
include { GUNZIP_REFFASTA   as GUNZIP_QUERYFASTA   } from '../../modules/local/gunzip_reffasta/main'
include { SAMTOOLS_FAIDX                           } from '../../modules/nf-core/samtools/faidx/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/gatk4/createsequencedictionary/main'


workflow PREPARE_REFERENCES {

    take:
    queryurl
    hosturl

    main:
    ch_queryref = channel.value(queryurl)
    ch_hostref  = channel.value(hosturl)


    ch_query    = DOWNLOAD_QUERYFASTA(ch_queryref).fasta
    if (hosturl) {
        hostfasta = DOWNLOAD_HOSTFASTA(ch_hostref).fasta
    } else {
        hostfasta = channel.value([])
    }
    queryfasta  = GUNZIP_QUERYFASTA(ch_query).gunzipped

    // These resources describe one reference, shared by every sample and interval.
    queryfasta  = queryfasta.map { ref -> tuple([id: ref.baseName], ref) }

    ch_faidx_in = queryfasta
                    .map { meta, fasta ->
                        return tuple(meta, fasta, [])
    }
    queryfai = SAMTOOLS_FAIDX(ch_faidx_in, false).fai
    querydict = GATK4_CREATESEQUENCEDICTIONARY(queryfasta).dict

    emit:
    queryfasta
    queryfai
    querydict
    hostfasta
}
