#!/bin/bash
#DATMA 
#Installation script
#Please see the DATMA manual for details
#

#Installing dependencies
#If you are using a different Linux distribution than Ubuntu
#please modify the next lines according your distribution
echo 'Installing dependencies'
sudo yum update -y
sudo yum install libz-dev
sudo yum install libboost-iostreams-dev
sudo yum install build-essential cmake
sudo yum install curl
sudo yum install ant
sudo yum install zlib1g-dev
sudo yum install -y pkg-config libfreetype6-dev libpng-dev python-matplotlib
sudo yum install numpy

#install BWA
sudo yum install -y bwa

#install SAM Tools
sudo yum install samtools

#install CheckM
sudo pip install checkm-genome
sudo pip install checkm-genome --upgrade --no-deps
sudo checkm data update

#Making the bin directory
cd tools
mkdir bin

#Trimmomatic
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip
unzip Trimmomatic-0.38.zip

#SDLS-lite 
echo 'Installing FM-index library'
git clone https://github.com/simongog/sdsl-lite.git
cd sdsl-lite
./install.sh
cd ..

#install selectFasta
echo 'Installing selectFasta'
cd selectFasta
make
cp selectFasta ../bin/
cd ..

#rdp-classifier
echo 'Installing rdp-classifier'
git clone https://github.com/rdpstaff/RDPTools.git
cd RDPTools/
git submodule init
git submodule update
make
cd ..

#install RAPIFILT
echo 'Installing RAPIFILT'
cd rapifilt
make
cp rapifilt ../bin/
cd ..

#install mergeNotCombined
echo 'Installing mergeNotCombined'
cd mergeNotCombined
make
cp mergeNotCombined ../bin/
cd ..

#install mapping (SDSL installed)
echo 'Installing mapping'
cd mapping
make
cp mapping ../bin/
cd ..

#genFm9
echo 'Installing genFM9'
cd genFm9
make
cp genFm9 ../bin/
cd ..

#Flash
echo 'Installing Flash'
git clone https://github.com/dstreett/FLASH2.git
cd FLASH2
make
cp flash2 ../bin/flash
cd ..

#CLAME
echo 'Installing CLAME'
cd CLAME
make
cp clame ../bin/
cd ..

#Megahit
echo 'Installing Megahit'
git clone https://github.com/voutcn/megahit.git
cd megahit
make
cp megahit ../bin/ #posiblemente toque pegar los otros binarios
cd ..

#SPAdes
echo 'Installing SPAdes'
wget http://cab.spbu.ru/files/release3.13.0/SPAdes-3.13.0.tar.gz
tar -xzf SPAdes-3.13.0.tar.gz
cd SPAdes-3.13.0
sudo ./spades_compile.sh
pwd=`pwd`
ln -s $pwd/spades.py ../bin/
cd ..

#Velvet
echo 'Installing Velvet'
git clone https://github.com/dzerbino/velvet.git
cd velvet
make
cp velvet* ../bin/
cd ..

#QUAST
echo 'Installing QUAST tool'
wget https://downloads.sourceforge.net/project/quast/quast-5.0.2.tar.gz
tar -xzf quast-5.0.2.tar.gz
cd quast-5.0.2
sudo ./setup.py install_full
cd ..

#Prodigal
echo 'Installing Prodigal'
git clone https://github.com/hyattpd/Prodigal.git
cd Prodigal/
make
cp prodigal ../bin/
cd ..

#krona 
echo 'Installing Krona'
wget https://github.com/marbl/Krona/releases/download/v2.7/KronaTools-2.7.tar
tar -xvf  KronaTools-2.7.tar
cd KronaTools-2.7/
sudo ./install.pl
./updateTaxonomy.sh
./updateAccessions.sh
cd ..

#Kaiju
echo 'Installing Kaiju'
git clone https://github.com/bioinformatics-centre/kaiju.git
cd kaiju/src/
make
cd ..
read -p "Do you want to download the Kaiju database (~27GB)?: (Y/N) " ans   
if [ "$ans" = 'Y' ] || [ "$ans" = 'y' ]; then
    echo "Downloading locat kaijudb"
    mkdir kaijudb
    cd kaijudb
    ../bin/makeDB.sh -r
    rm -r genomes/
    rm kaiju_db.bwt
    rm kaiju_db.faa
    rm kaiju_db.sa
    cd ..
else
    echo "You can specify the Kaiju database using the configuration file"
fi
cp bin/* ../bin/
cd ..

#Blast
echo 'Installing BLAST'
sudo yum install ncbi-blast+

#update nt database
read -p "Do you want to download the NT-database (~58GB)?: (Y/N) " ans   
if [ "$ans" = 'Y' ] || [ "$ans" = 'y' ]; then
    echo "Downloading locat NT"
    git clone https://github.com/jrherr/bioinformatics_scripts.git
    cp bioinformatics_scripts/perl_scripts/update_blastdb.pl bin/
    chmod +x bin/update_blastdb.pl
    cd ../
    mkdir blastdb
    cd blastdb
    perl ../tools/bin/update_blastdb.pl --passive --decompress nt
    cd ..
else
    echo "You can specify a local NT-database using the configuration file"
fi

#FINISH
echo 'DATMA INSTALLED'
echo 'Please update your PATH to include build tools'
echo 'export PATH=<intall_PATH>/DATMA/tools/bin/:$PATH'
