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
raw_files = file("/home/robbe/ionbot/${params.accession}/raw_files/")
raw_files_dir = raw_files.mkdirs()
println raw_files_dir ? "Directory $raw_files created" : "Could not create directory: $raw_files_dir"
mgf_files = file("/home/robbe/ionbot/${params.accession}/mgf_files")
mgf_files_dir = mgf_files.mkdirs()
println mgf_files_dir ? "Directory $mgf_files created" : "Could not create directory: $mgf_files_dir"
ionbot_files = file("/home/robbe/ionbot/${params.accession}/ionbot_files")
ionbot_files_dir = ionbot_files.mkdirs()
println ionbot_files_dir ? "Directory $ionbot_files created" : "Could not create directory: $raw_files_dir"
configfile = file("/home/robbe/ionbot/configfiles/config${params.db}.txt")
println "Config file used: $configfile"

//Channels
Rawchannel = Channel.fromPath("/home/robbe/ionbot/${params.accession}/raw_files/*.raw")
MGFchannel = Channel.fromPath("/home/robbe/ionbot/${params.accession}/mgf_files/*.mgf")

//Download raw files from PRIDE
process PRIDE_download {

    """
    pridepy.py download-all-raw-files -a ${params.accession} -o /home/robbe/ionbot/${params.accession}/raw_files/
    """  
}


//Converts raw files to mgf files
process raw_to_mgf {

    container 'quay.io/biocontainers/thermorawfileparser:1.2.3--1'

    input:
    file rawFile from Rawchannel

    """
    ThermoRawFileParser -i=${rawFile} -m=0 -f=0 -o=/home/robbe/ionbot/${params.accession}/mgf_files --ignoreInstrumentErrors
    """
}



//search mgf files with ionbot
process ionbotsearch {

    container 'gcr.io/omega-cloud-195908/ionbot:v0.8.0'
    cpus 16

    input:
    file mgfFile from MGFchannel

    """
    ionbot -c "${configfile}" -o "/home/robbe/ionbot/${params.accession}/ionbot_files" -a ${task.cpus} -I -R -e -m -r -X "${mgfFile}"
    echo "The MGF file ${mgfFile} has been succesfully searched with ionbot"
    """
}

//Removing the hash from the file names caused by pridepy
Files_to_rename = file("/home/robbe/ionbot/${params.accession}/ionbot_files")
until_hyphen = ~/.*(?=-)/
filestorename = Files_to_rename.listFiles()
for( def file : filestorename ) {
    name = file.getName()
    new_name = name - until_hyphen - "-"
    file.renameTo("/home/robbe/ionbot/${params.accession}/ionbot_files/$new_name")}
