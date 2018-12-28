---------------------------------------------------------------
CLAME:"Clasificador Metagenomico"
---------------------------------------------------------------
Description:
    CLAME is a binning software for metagenomic reads.
    It immplements a fm-index search algorithm for nucleotide 
    sequence alignment. Then it uses strongly connected component strategy
    to bin sequences with similar DNA composition.

---------------------------------------------------------------
Installation
---------------------------------------------------------------
1. To download and install the SDSL - Succinct Data Structure Library, follow the instructions on the web page.
https://github.com/simongog/sdsl-lite

2. To download CLAME source codes
git clone https://github.com/andvides/CLAME.git

3. cd CLAME/

4. make

5. Type ./clame -h to show the help

---------------------------------------------------------------
CLAME Versions
---------------------------------------------------------------
1. Type make or make all: To Compile CLAME.
2. Type make debug: To compile a debug version of CLAME.

---------------------------------------------------------------
Usage
---------------------------------------------------------------
./clame -b 70 -multiFasta test/Bancomini.fna -nt 4 -output bminiBins -print

Output files
1. bminiBins_*.fasta:   Output fasta file for all the bins reported
2. bminiBins.binning:   All bins reported 
3. bmini.fm9:           FM-index output
4. bmini.index:         First colum contains the origal name for each read, the second column the index used by CLAME
5. bminiBins.links:     Histogram links by number of reads
6. bminiBins.result:    Adjacency list for the overlap detected by each read
     
---------------------------------------------------------------
Authors
---------------------------------------------------------------
Benavides A(1), Alzate JF (2),(3) and Cabarcas F (1),(3)
1.	Grupo Sistemic, Departamento de Ingeniería Electrónica, Facultad de Ingenieria, Universidad de Antioquia.
2.	Centro Nacional de Secuenciacion Genomica-CNSG, Sede de Investigación Universitaria SIU, Universidad de Antioquia
3.	Grupo de Parasitología, Departamento de Microbiología y Parasitología, Facultad de Medicina, Universidad de Antioquia
