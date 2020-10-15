#!/bin/sh
### install MA57
### Note: if using gcc compilers, running this script may give errors
### like "Reference to unimplemented intrinsic" if g77 is installed.
### Run using F77=gfortran ./installMa57.sh instead.

#assume ma57 in tar.gz file
fn=`ls ma57*.tar.gz`
name=`basename ${fn} .tar.gz`
tar -zxvf $fn
ln -s ./${name} ./src

FULLPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

cd src
./configure FFLAGS=-fPIC --with-metis=${FULLPATH}/../METIS_4/src/libmetis.a --prefix=`pwd`

make clean
make install
make check
