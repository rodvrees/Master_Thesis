#!/usr/bin/env nextflow

// PRIDE accession number
params.accession = "PXD014381"

rawpath = "/public/compomics2/Robbe/${params.accession}/*.raw"

Rawchannel = Channel.fromPath("$rawpath")
//Convert raw files to mgf files
process raw_to_mgf {

    container 'quay.io/biocontainers/thermorawfileparser:1.2.3--1'

    input:
    each path(rawFile) from Rawchannel

    """
    ThermoRawFileParser -i=${rawFile} -m=0 -f=1 -o="/public/compomics2/Robbe/${params.accession}/" --ignoreInstrumentErrors
    """
}
