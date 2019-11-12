#!/bin/bash
#DATMA 
#Installation script
#Please see the DATMA manual for details
#

#Making the bin directory
cd tools
mkdir bin

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

#genFm9
echo 'Installing binning'
cd binning
make
cp binning ../bin/
cd ..

