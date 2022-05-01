#!/usr/bin/env nextflow

accession = Channel.from("PXD014629")

process PRIDE_download {

    input: 
    val x from accession

    """
    pridepy.py download-all-raw-files -a ${x} -o /home/robbe/ionbot/${x}/raw_files/
    """  
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
