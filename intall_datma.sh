#Download the DATMA github
git clone https://github.com/dstreett/DATMA.git
./install_datma.sh

#Install dependency
sudo apt-get install libz-dev
sudo apt-get install curl

#SDLS-lite 
git clone https://github.com/simongog/sdsl-lite.git
cd sdsl-lite
./install.sh
cd ..

#Making the bin directory
cd tools
mkdir bin

#install selectFasta
cd selectFasta
make
cp selectFasta ../bin/
cd ..

#install RAPIFILT
cd tools
cd rapifilt
make
cp rapifilt ../bin/
cd ..


#install mapping (SDSL installed)
cd tools
cd mapping
make
cp mapping ../bin/
cd ..

#genFm9
cd tools
cd genFm9
make
cp genFm9 ../bin/
cd ..

#Flash
git clone https://github.com/dstreett/FLASH2.git
cd FLASH2
make
cp flash2 ../bin/flash

#CLAME
cd tools
cd CLAME
make
cp clame ../bin/
cd ..


#Megahit
git clone https://github.com/voutcn/megahit.git
cd megahit
make
cp megahit ../bin/ #posiblemente toque pegar los otros binarios

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
cp bin/* ../bin/
cd ..

#Prodigal
https://github.com/hyattpd/Prodigal.git
cd Prodigal/
make
cp prodigal ../bin/

#krona 
wget https://github.com/marbl/Krona/releases/download/v2.7/KronaTools-2.7.tar
tar -xvf  KronaTools-2.7.tar
cd KronaTools-2.7/
sudo ./install.pl
./updateTaxonomy.sh
./updateAccessions.sh



