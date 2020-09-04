#!/bin/bash
if [ ! -f "metis-4.0.3.tar.gz" ]; then
  echo " "
  echo "##### Downloading METIS 4"
  echo " "

  fn=metis-4.0.3.tar.gz
  echo "### Downloading Metis:"
  if wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/${fn}
  then
    echo "### Metis4: Download Successful.\n"
  else
    echo "### Metis4: Download Failed.\n"
    exit 1
  fi
fi

if [ ! -d "./metis-4.0.3" ]; then
  name=`basename ${fn} .tar.gz`
  tar -zxvf $fn
  ln -s ./${name} ./src
fi

#compile metis
echo "### Compiling METIS 4"

cp ./Makefile.in src/Makefile.in

cd src
make clean
make

cd Graphs
./mtest ./test.mgraph
./mtest ./4elt.graph

cd ../..
