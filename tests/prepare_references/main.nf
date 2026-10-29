// Isolated harness: no pipeline plugins, samplesheet or downstream analysis.
include { PREPARE_REFERENCES } from '../../subworkflows/local/prepare_references'

params.query_fasta = null
params.host_fasta = null
params.hosturl = null

workflow {
    PREPARE_REFERENCES(params.queryurl, params.hosturl, params.query_fasta, params.host_fasta)
    // Singleton queues would silently drop all but one sample downstream.
    assert PREPARE_REFERENCES.out.queryfasta instanceof groovyx.gpars.dataflow.DataflowVariable
    assert PREPARE_REFERENCES.out.queryfai instanceof groovyx.gpars.dataflow.DataflowVariable
    assert PREPARE_REFERENCES.out.querydict instanceof groovyx.gpars.dataflow.DataflowVariable
    assert PREPARE_REFERENCES.out.hostfasta instanceof groovyx.gpars.dataflow.DataflowVariable
    PREPARE_REFERENCES.out.queryfasta.view { meta, fasta ->
        "QUERY\t${meta.id}\t${fasta.name}\t${fasta.text.trim()}"
    }
    PREPARE_REFERENCES.out.queryfai.view { meta, fai ->
        "FAI\t${meta.id}\t${fai.name}\t${fai.text.trim()}"
    }
    PREPARE_REFERENCES.out.querydict.view { meta, dict ->
        assert dict.text.contains('@SQ\tSN:synthetic\tLN:')
        "DICT\t${meta.id}\t${dict.name}"
    }
    PREPARE_REFERENCES.out.hostfasta.view { fasta ->
        fasta ? "HOST\t${fasta.name}\t${fasta.text.trim()}" : 'HOST\tabsent'
    }
}
