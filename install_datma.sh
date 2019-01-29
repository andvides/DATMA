#!/bin/bash
#DATMA 
#Installation script
#Please see the DATMA manual for details
#

#Install dependency
#If you are using a different Linux distribution than Ubuntu
#please modify the next lines according your distribution
sudo apt-get install libz-dev
sudo apt-get install libboost-iostreams-dev
sudo apt-get install build-essential cmake
sudo apt-get install curl
sudo apt install ant

#Making the bin directory
cd tools
mkdir bin


#SDLS-lite 
git clone https://github.com/simongog/sdsl-lite.git
cd sdsl-lite
./install.sh
cd ..

#install selectFasta
cd selectFasta
make
cp selectFasta ../bin/
cd ..

#rdp classifier
git clone https://github.com/rdpstaff/RDPTools.git
cd RDPTools/
git submodule init
git submodule update
make
cd ..

#install RAPIFILT
cd rapifilt
make
cp rapifilt ../bin/
cd ..


#install mapping (SDSL installed)
cd mapping
make
cp mapping ../bin/
cd ..

#genFm9
cd genFm9
make
cp genFm9 ../bin/
cd ..

#Flash
git clone https://github.com/dstreett/FLASH2.git
cd FLASH2
make
cp flash2 ../bin/flash
cd ..

#CLAME
cd CLAME
make
cp clame ../bin/
cd ..

#Megahit
git clone https://github.com/voutcn/megahit.git
cd megahit
make
cp megahit ../bin/ #posiblemente toque pegar los otros binarios
cd ..

#SPAdes
wget http://cab.spbu.ru/files/release3.13.0/SPAdes-3.13.0.tar.gz
tar -xzf SPAdes-3.13.0.tar.gz
cd SPAdes-3.13.0
sudo ./spades_compile.sh
cd ..

#Velvet
git clone https://github.com/dzerbino/velvet.git
cd velvet
make
cp velvet* ../bin/
cd ..

#Blast
sudo apt-get install ncbi-blast+

#Kaiju
git clone https://github.com/bioinformatics-centre/kaiju.git
cd kaiju/src/
make
cd ..
mkdir kaijudb
cd kaijudb
../bin/makeDB.sh -r
rm -r genomes/
rm kaiju_db.bwt
rm kaiju_db.faa
rm kaiju_db.sa
cd ..
cp bin/* ../bin/
cd ..

#Prodigal
git clone https://github.com/hyattpd/Prodigal.git
cd Prodigal/
make
cp prodigal ../bin/
cd ..

#krona 
wget https://github.com/marbl/Krona/releases/download/v2.7/KronaTools-2.7.tar
tar -xvf  KronaTools-2.7.tar
cd KronaTools-2.7/
sudo ./install.pl
./updateTaxonomy.sh
./updateAccessions.sh
cd ..

#update nt database
git clone https://github.com/jrherr/bioinformatics_scripts.git
cp bioinformatics_scripts/perl_scripts/update_blastdb.pl bin/
chmod +x bin/update_blastdb.pl
cd ../
mkdir blastdb
cd blastdb
perl ../tools/bin/update_blastdb.pl --passive --decompress nt
cd ..
#FINISH
echo 'DATMA INSTALLED'
