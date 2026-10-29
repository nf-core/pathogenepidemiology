include { PREPARE_REFERENCES } from '../main.nf'

process REUSE_REFERENCE {
    input:
    val sample
    tuple val(meta), path(fasta)
    tuple val(fai_meta), path(fai)
    tuple val(dict_meta), path(dict)
    path host

    output:
    tuple val(sample), path('checked.txt'), emit: checked

    script:
    """
    test -s ${fasta}
    test -f ${fai}
    test -f ${dict}
    echo '${sample}:${host ? "host" : "no-host"}' > checked.txt
    """
}

workflow REFERENCE_VALUES {
    take:
    queryurl
    hosturl

    main:
    PREPARE_REFERENCES(queryurl, hosturl, null, null)
    // Assert channel type, not just a singleton emission: queues also emit once.
    assert PREPARE_REFERENCES.out.queryfasta instanceof groovyx.gpars.dataflow.DataflowVariable
    assert PREPARE_REFERENCES.out.queryfai instanceof groovyx.gpars.dataflow.DataflowVariable
    assert PREPARE_REFERENCES.out.querydict instanceof groovyx.gpars.dataflow.DataflowVariable
    assert PREPARE_REFERENCES.out.hostfasta instanceof groovyx.gpars.dataflow.DataflowVariable
    REUSE_REFERENCE(channel.of('sample1', 'sample2'),
        PREPARE_REFERENCES.out.queryfasta, PREPARE_REFERENCES.out.queryfai,
        PREPARE_REFERENCES.out.querydict, PREPARE_REFERENCES.out.hostfasta)

    emit:
    checked = REUSE_REFERENCE.out.checked
    queryfasta = PREPARE_REFERENCES.out.queryfasta
    queryfai = PREPARE_REFERENCES.out.queryfai
    querydict = PREPARE_REFERENCES.out.querydict
    hostfasta = PREPARE_REFERENCES.out.hostfasta
}
