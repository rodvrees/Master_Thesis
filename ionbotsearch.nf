#!/usr/bin/env nextflow

// PRIDE accession number
params.accession = "PXD014629"

/*
Database: Should be either human, mouse, humanTMT10, mouseTMT10 (and other, make configfiles accordingly), 
alternatively, make configfiles and make the parameter the path, maybe even better
default = human
*/
params.db = "human"

//Create needed directories
ionbot_files = file("/home/robbe/ionbot/${params.accession}/ionbot_files")
ionbot_files_dir = ionbot_files.mkdirs()
println ionbot_files_dir ? "Directory $ionbot_files created" : "Could not create directory: $raw_files_dir"
configfile = file("/home/robbe/ionbot/configfiles/config${params.db}.txt")
println "Config file used: $configfile"

//Download raw files from PRIDE
process PRIDE_download {
    output:
    path '*.raw' into Rawchannel
    """
    pridepy.py download-all-raw-files -a ${params.accession} -o $task.workDir
    """  
}

Rawchannel.subscribe onNext: { println "File: ${it.name}" }

//Converts raw files to mgf files
process raw_to_mgf {

    container 'quay.io/biocontainers/thermorawfileparser:1.3.4--ha8f3691_0'

    input:
    each path(rawFile) from Rawchannel
    output:
    path '*.mgf' into MGFchannel

    """
    ThermoRawFileParser -i=${rawFile} -m=0 -f=0 --ignoreInstrumentErrors
    """
}


//search mgf files with ionbot
process ionbotsearch {

    container 'gcr.io/omega-cloud-195908/ionbot:v0.9.0'
    cpus 16

    input:
    each mgfFile from MGFchannel
    file(configfile)

    """
    mkdir /home/robbe/ionbot/${params.accession}/ionbot_files/${mgfFile.baseName}
    ionbot -c "${configfile}" -o "/home/robbe/ionbot/${params.accession}/ionbot_files/${mgfFile.baseName}" -a ${task.cpus} -I -R -e -m -r "${mgfFile}"
    """
}


    