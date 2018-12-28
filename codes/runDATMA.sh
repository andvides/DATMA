#rm -r ../examples/crypto/bins/
#mv ../examples/bmini/round_0_b100/readsForbin.fastq examples/crypto/
#rm -r ../examples/crypto/round_*
#rm -r ../examples/bmini
#rm -r ~/DATMA_*
#COMPSs enviroment
#export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64/
#export PATH=$PATH:$JAVA_HOME/bin
#export CLASSPATH=$JAVA_HOME/jre/lib/ext:$JAVA_HOME/lib/tools.jar

#runcompss --master_name=pc -d --lang=python /home/andrespc/Documents/datma/codes/datma.py -f $1
#python /home/users/andresb/datma/codes/finalReport.py -f $1

python datma_seq.py -f ../examples/configBmini.txt
