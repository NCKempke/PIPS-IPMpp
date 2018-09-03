#!/bin/bash

# set default values
gams_file="eps"
nblocks="3"
np="1"
scale=""
stepLp=""
presolve=""

for i in "$@"
do
case $i in
    -GAMSFILE=*|--GAMSFILE=*)
    gams_file="${i#*=}"
    shift # past argument=value
    ;;
    -BLOCKS=*|--BLOCKS=*)
    nblocks="${i#*=}"
    shift # past argument=value
    ;;
    -NP=*|--NP=*)
    np="${i#*=}"
    shift # past argument=value
    ;;
    -SCALE=*|--SCALE=*)
    scale="${i#*=}"
    shift # past argument=value
    ;;
    -STEPLP=*|--STEPLP=*)
    stepLp="${i#*=}"
    shift # past argument=value
    ;;
    -PRESOLVE=*|--PRESOLVE=*)
    presolve="${i#*=}"
    shift # past argument=value
    ;;
    *)
          # unknown option
    ;;
esac
done
    
if [ -f $gams_file.gms ]; then
   echo "File $FILE exists."
else
   echo "File $gams_file.gms does not exist. Please give a valid .gms file without the .gms extention. Exiting script."
   exit
fi

gams $gams_file --METHOD=PIPS > /dev/null
../../../../build_pips/gmschk -g $GAMSSYSDIR -T -X $nblocks $gams_file.gdx > /dev/null

if [ "$stepLp" = "true" ]; then
  stepLp="stepLp"
else
  stepLp=""
fi

if [ "$presolve" = "true" ]; then
  presolve="presolve"
else
  presolve=""
fi

if [ "$scale" = "true" ]; then
  scale="scale"
else
  scale=""
fi

mpirun -np $np ../../../../build_pips/gmspips $nblocks $gams_file $GAMSSYSDIR $scale $stepLp $presolve 2>&1 | tee pips.out


