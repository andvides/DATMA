# DATMA
Distributed AuTomatic Metagenomic Assembly and Annotation framework
--------------------------------------------------------------
Introduction
--------------------------------------------------------------
DATMA is a distributed automatic pipeline for fast metagenomic analysis that includes: sequencing quality control, 16S-identification, reads binning, de novo assembly, ORF detection and taxonomic annotation.

--------------------------------------------------------------
Installing DATMA
--------------------------------------------------------------
NON distributed version installation (NON COMPSs support)<br/>
1. Clone the DATMA directory:<br/> 
git clone https://github.com/andvides/DATMA.git
2. Run the ./install_datma.sh script:<br/>
bash install_datma.sh 
3. Export the PATH:<br/> 
export PATH=<install_path>/DATMA/tools/bin/:$PATH

Distributed version installation (COMPSs support)<br/>
1. Intall COMPSs framework (see the COMPSs manual)
2. Clone the DATMA directory:<br/> 
git clone https://github.com/andvides/DATMA.git
3. Run the ./install_datma.sh script:<br/>
bash install_datma.sh 
4. Export the PATH:<br/>
export PATH=<install_path>/DATMA/tools/bin/:$PATH

---------------------------------------------------------------
Running DATMA
--------------------------------------------------------------
1. Generate the 16S database index:<br/>
see <install_path>/DATMA/16sDatabases/README
2. Edit the configBmini.txt:<br/>
nano <install_path>/DATMA/examples/configBmini.txt
3. Run datma and select version (seq) or (compss):<br/>
<install_path>/DATMA/runDATMA.sh <install_path>/datma/examples/configBmini.txt seq

---------------------------------------------------------------
Output
---------------------------------------------------------------
1. 16sSeq: Directory with the identified 16S sequences
2. bins: Assembly and annotation of the bins
3. clean: Quality control report
4. round: CLAME executions
5. readsForbin.fastq: Balance of reads without classification
6. binsBlastn.html: Blastn annotation for the contigs in Krona format
7. binsKaiju.html: Kaiju annotation for the contigs in Krona format
8. report.html: Bins report
9. resumenFile.html: Darma summarized report

---------------------------------------------------------------
Authors
---------------------------------------------------------------
Benavides A(1), Sanchez F (2), Alzate JF (3),(4) and Cabarcas F (1),(3)

1. Grupo Sistemic, Departamento de Ingeniería Electrónica, Facultad de Ingenieria, Universidad de Antioquia.
2. Smart Variable S.L, Barcelona, Spain
3. Centro Nacional de Secuenciacion Genomica-CNSG, Sede de Investigación Universitaria SIU, Universidad de Antioquia
4. Grupo de Parasitología, Departamento de Microbiología y Parasitología, Facultad de Medicina, Universidad de Antioquia

---------------------------------------------------------------
FAQ
--------------------------------------------------------------
Please contact us:<br/>
bernardo.benavides@udea.edu.co<br/>
felipe.cabarcas@udea.edu.co
