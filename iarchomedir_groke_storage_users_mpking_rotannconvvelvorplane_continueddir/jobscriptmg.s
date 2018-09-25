#!/bin/bash
#PBS -q standard         
#PBS -l select=1:ncpus=4:node_type=4way
#PBS -l walltime=36:00:00 
#PBS -j oe
               
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=4
export PSC_OMP_AFFINITY=TRUE

/usr/bin/time ./velvor2dnewsch3.4.exe < ./input8e78e4pr1-1r1.s
