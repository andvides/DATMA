#!/bin/bash
chmod +x datma

mkdir -p $PREFIX/bin/
cp codes/*py $PREFIX/bin/
cp datma $PREFIX/bin/
echo "export datmaPATH=$PREFIX/bin/" >> ~/.basrc
