#!/bin/sh
### install MC30

#we assume MC30 is in a tar.gz file
fn=`ls mc30*.tar.gz`
name=`basename ${fn} .tar.gz`
tar -zxvf $fn
ln -s ./${name} ./src

#configure and build mc30
cd src
./configure FFLAGS=-fPIC --prefix=`pwd`
make -j4 install
make check
