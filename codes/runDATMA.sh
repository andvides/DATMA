#!/bin/bash
#COMPSs enviroment
#export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64/ #PATH to JAVA HOME
#export PATH=$PATH:$JAVA_HOME/bin
#export CLASSPATH=$JAVA_HOME/jre/lib/ext:$JAVA_HOME/lib/tools.jar

config="$1" #configuration file
mode="$2" #mode sequential (seq) or distributed (compss)

if [ "$2" = 'seq' ]; then
    echo "Running DATMA in sequential mode."
    python datma_seq.py -f $config
elif [ "$2" = 'compss' ]; then
    echo "Running DATMA with COMPSs support."
    runcompss --master_name=pc -d --lang=python ./datma.py -f $config
else
    echo "Select between sequential or  distributed mode (see DATMA manual)"
    exit 1;
fi
