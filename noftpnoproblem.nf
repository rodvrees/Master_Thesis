#!/usr/bin/env nextflow

// PRIDE accession number
params.accession = "PXD014629"

/*
Database: Should be either human, mouse, humanTMT10, mouseTMT10 (and other, make configfiles accordingly), 
alternatively, make configfiles and make the parameter the path, maybe even better
default = human
*/
params.configpath = "$launchDir/configfiles/confighuman.txt"

//Create needed directories
ionbot_files = file("./${params.accession}/ionbot_files")
ionbot_files_dir = ionbot_files.mkdirs()
println ionbot_files_dir ? "Directory $ionbot_files created" : "Could not create directory: $ionbot_files_dir"
configfile = file("${params.configpath}")
println "Config file used: $configfile"

Rawchannel = Channel.fromPath("$launchDir/$params.accession/raw_files/*.raw")
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
    mkdir $launchDir/${params.accession}/ionbot_files/${mgfFile.baseName} -p
    ionbot -c "${configfile}" -o "$launchDir/${params.accession}/ionbot_files/${mgfFile.baseName}" -a ${task.cpus} -I -R -e -m -r "${mgfFile}"
    """
}


    