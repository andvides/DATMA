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
git clone https://github.com/dstreett/DATMA.git
2. running the intall script:
./install_datma.sh 
3. exporting the PATH: 
export PATH=<install_path>/datma/tools/bin/:$PATH

A full intallation can be done by (COMPSs support)
1. installing COMPSs framework (see the COMPSs manual)
2. clonning the DATMA directory 
git clone https://github.com/dstreett/DATMA.git
3. running the ./install_datma.sh script
./install_datma.sh 
4. exporting the PATH:
export PATH=<install_path>/datma/tools/bin/:$PATH

---------------------------------------------------------------
FAQ
---------------------------------------------------------------
Please contact us:
bernardo.benavides@udea.edu.co
felipe.cabarcas@udea.edu.co
