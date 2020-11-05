#!/bin/sh
### install MA60

#we assume MA60 is in a tar.gz file
fn=`ls ma60*.tar.gz`
name=`basename ${fn} .tar.gz`
tar -zxvf $fn
ln -s ./${name} ./src

cd src
./configure FFLAGS=-fPIC --prefix=`pwd`
make -j4 install
make check
