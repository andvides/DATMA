#!/bin/bash
#DATMA 
#Installation script
#Please see the DATMA manual for details
#

#Making the bin directory
cd tools
mkdir bin

#install selectFasta
echo 'Installing selectFasta'
cd selectFasta
make
cp selectFasta ../bin/
cd ..



