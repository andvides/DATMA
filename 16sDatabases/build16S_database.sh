#!/bin/bash
dataBase="$1" #path to the database: ~/DATMA/16sDatabases/rdp
mode="$2" #mode genFM9 (fm9) or BWA (bwa)

unzip "$dataBase.zip"
if [ "$2" = 'fm9' ]; then
    echo "Generating the FM-index"
    genFm9 -multiFasta "$dataBase.fna" -output $dataBase
elif [ "$2" = 'bwa' ]; then
    echo "Generating the BWA index"
    bwa index  "$dataBase.fna"
else
    echo "Select between FM-index or BWA index (see DATMA manual)"
    exit 1;
fi
