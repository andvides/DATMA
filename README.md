# DATMA
Distributed AuTomatic Metagenomic Assembly and Annotation framework
---------------------------------------------------------------
Introduction
---------------------------------------------------------------
DATMA is a distributed automatic pipeline for fast metagenomic analysis that includes: sequencing quality control, 16S-identification, reads binning, de novo assembly, ORF detection and taxonomic annotation.

---------------------------------------------------------------
Install
---------------------------------------------------------------
A bassic intallation can be done by (NO COMPSs support)
1. clonning the DATMA directory: 
git clone https://github.com/andvides/DATMA.git
2. running the intall script:
./install_datma.sh 
3. exporting the PATH: 
export PATH=<install_path>/datma/tools/bin/:$PATH

A full intallation can be done by (COMPSs support)
1. installing COMPSs framework (see the COMPSs manual)
2. clonning the DATMA directory: 
git clone https://github.com/andvides/DATMA.git
3. running the ./install_datma.sh script:
./install_datma.sh 
4. exporting the PATH:
export PATH=<install_path>/datma/tools/bin/:$PATH

---------------------------------------------------------------
Running
---------------------------------------------------------------
1. Generate the 16S database index
cat <install_path>/datma/16sDatabases/README
2. Edit the configBmini.txt
nano <install_path>/datma/examples/configBmini.txt
3. Run datma
<install_path>/datma/runDATMA.sh <install_path>/datma/examples/configBmini.txt

---------------------------------------------------------------
Output
---------------------------------------------------------------

16sSeq: Directory with the identified 16S sequences
bins: Assembly and annotation of the bins
clean: Quality control report
round: CLAME executions
readsForbin.fastq: Balance of reads without classification
binsBlastn.html: Blastn annotation for the contigs in Krona format
binsKaiju.html: Kaiju annotation for the contigs in Krona format
report.html: Bins report
resumenFile.html: Darma summarized report

---------------------------------------------------------------
FAQ
---------------------------------------------------------------
# DATMA
Distributed AuTomatic Metagenomic Assembly and Annotation framework
---------------------------------------------------------------
Introduction
---------------------------------------------------------------
DATMA is a distributed automatic pipeline for fast metagenomic analysis that includes: sequencing quality control, 16S-identification, reads binning, de novo assembly, ORF detection and taxonomic annotation.

---------------------------------------------------------------
Install
---------------------------------------------------------------
A bassic intallation can be done by (NO COMPSs support)
1. clonning the DATMA directory: 
git clone https://github.com/andvides/DATMA.git
2. running the intall script:
./install_datma.sh 
3. exporting the PATH: 
export PATH=<install_path>/datma/tools/bin/:$PATH

A full intallation can be done by (COMPSs support)
1. installing COMPSs framework (see the COMPSs manual)
2. clonning the DATMA directory: 
git clone https://github.com/andvides/DATMA.git
3. running the ./install_datma.sh script:
./install_datma.sh 
4. exporting the PATH:
export PATH=<install_path>/datma/tools/bin/:$PATH

---------------------------------------------------------------
Running
---------------------------------------------------------------
1. Generate the 16S database index
cat <install_path>/datma/16sDatabases/README
2. Edit the configBmini.txt
nano <install_path>/datma/examples/configBmini.txt
3. Run datma
<install_path>/datma/runDATMA.sh <install_path>/datma/examples/configBmini.txt

---------------------------------------------------------------
Output
---------------------------------------------------------------
16sSeq: Directory with the identified 16S sequences
bins: Assembly and annotation of the bins
clean: Quality control report
round: CLAME executions
readsForbin.fastq: Balance of reads without classification
binsBlastn.html: Blastn annotation for the contigs in Krona format
binsKaiju.html: Kaiju annotation for the contigs in Krona format
report.html: Bins report
resumenFile.html: Darma summarized report

---------------------------------------------------------------
FAQ
---------------------------------------------------------------
Please contact us:
bernardo.benavides@udea.edu.co
felipe.cabarcas@udea.edu.co

