#!/bin/sh

if [ -d "$HOME/TTBarAnalysis/build" ]; then
   rm -rf "$HOME/TTBarAnalysis/build"
fi

mkdir "$HOME/TTBarAnalysis/build"
cd "$HOME/TTBarAnalysis/build"
cmake -C $ILCSOFT/ILCSoft.cmake ..
make
make install
