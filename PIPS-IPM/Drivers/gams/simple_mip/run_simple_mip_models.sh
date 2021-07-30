#!/bin/bash
# unpack, split and run simple_mip_models

echo "Running Simple-MIP"

n_procs=1
if [ "$#" -ne 1 ]; then
    echo "Expecting excatly one input argument specifying the number of mpi-processes to run - none given!"
    echo "Defaulting to #procs = 1"
else
  n_procs=$1
fi

mpirun_command=mpirun
mpi_nprocs_flag=np
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

if [ -d "${DIR}/5_blk" ]
then
    echo "Removing old blockfiles"
    rm -r ${DIR}/*_blk
fi

echo "Unpacking Archive..."
unzip $DIR/mip_gdx_files.zip -d $DIR/

echo "Splitting models using gmschk"
mkdir $DIR/5_blk
mkdir $DIR/56_blk
mkdir $DIR/366_blk

mv $DIR/allblocks_mip_5blk* $DIR/5_blk
mv $DIR/allblocks_mip_56blk* $DIR/56_blk
mv $DIR/allblocks_mip_366blk* $DIR/366_blk

$DIR/../../../../build/gmschk -g $GAMSSYSDIR -T -X 5 $DIR/5_blk/allblocks_mip_5blk.gdx
$DIR/../../../../build/gmschk -g $GAMSSYSDIR -T -X 56 $DIR/56_blk/allblocks_mip_56blk.gdx
$DIR/../../../../build/gmschk -g $GAMSSYSDIR -T -X 366 $DIR/366_blk/allblocks_mip_366blk.gdx

echo "Running the models one by one"

min4=$(( ${n_procs} <= 4 ? ${n_procs} : 4 ))
echo "Running ${mpirun_command} -${mpi_nprocs_flag} ${min4} ${DIR}/../../../../build/gmspips 5 ${DIR}/5_blk/allblocks_mip_5blk ${GAMSSYSDIR} presolve scaleGeo stepLp"
$mpirun_command -$mpi_nprocs_flag ${min4} $DIR/../../../../build/gmspips 5 $DIR/5_blk/allblocks_mip_5blk $GAMSSYSDIR presolve scaleGeo stepLp

min55=$(( ${n_procs} <= 55 ? ${n_procs} : 55 ))
echo "Running ${mpirun_command} -${mpi_nprocs_flag} ${min55} ${DIR}/../../../../build/gmspips 56 ${DIR}/56_blk/allblocks_mip_56blk ${GAMSSYSDIR} presolve scaleGeo stepLp"
$mpirun_command -$mpi_nprocs_flag ${min55} $DIR/../../../../build/gmspips 56 $DIR/56_blk/allblocks_mip_56blk $GAMSSYSDIR presolve scaleGeo stepLp

min365= $(( ${n_procs} < 366 ? ${n_procs} : 365 ))
echo "Running ${mpirun_command} -${mpi_nprocs_flag} ${min365} ${DIR}/../../../../build/gmspips 366 ${DIR}/366_blk/allblocks_mip_366blk ${GAMSSYSDIR} presolve scaleGeo stepLp"
$mpirun_command -$mpi_nprocs_flag ${min365} $DIR/../../../../build/gmspips 366 $DIR/366_blk/allblocks_mip_366blk $GAMSSYSDIR presolve scaleGeo stepLp
