#!/bin/bash
#----------------------------------------------------------------------------------------------------------------------------------
# 					            SRA Files Downloading
#----------------------------------------------------------------------------------------------------------------------------------
# This script takes as input SRA files from GEO experiments and obtains the corresponding FASTQ files.
# It must have been installed the SRA Toolkit.

# Cheking and storage of the Arguments
USAGE="USAGE= $0 -d <date> -s <sra> | -s <sra,sra>"

if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
	echo "$USAGE"
	exit 1
fi
if [ "$#" -ne "4" ]; then
	echo "Incorrect arguments: $USAGE"
	exit 1
fi

while [[ $# -gt 0 ]]
do
key="$1"
case $key in
	-s|--sra)
	SRAs="$2"
	shift
	shift
	;;
	-d|--date)
	Date="$2"
	shift
	shift
	;;
esac
done

# Settings. ChIP-Seq identifiers: cell line and data
HomeDir=~/Bioinformatica_Núria/
mkdir $HomeDir/SRAs_$Date
FilesDir=$HomeDir/SRAs_$Date
SRAs=$(echo "$SRAs" | sed 's/,/ /')

#----------------------------------------------------------------------------------------------------------------------------------
# A for loop will iterate over all the files included in the working directory with the .gz extension (FASTQ files)
cd $FilesDir
for SRA in $SRAs
	do echo "-------------------- Processing sample: $SRA ---------------------"
	$HomeDir/sratoolkit.2.9.6-1-centos_linux64/bin/fastq-dump --split-spot $SRA
done

