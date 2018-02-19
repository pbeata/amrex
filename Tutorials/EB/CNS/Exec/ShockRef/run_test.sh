#!/bin/bash

inFile="./inputs.amr"

#outFile="./results/SF_4MPI_8OMP_10steps.txt"
#outFile="./results/KS_4MPI_8OMP_10steps.txt"
#outFile="./results/RR_4MPI_8OMP_10steps.txt"

export OMP_NUM_THREADS=8

make

#mpiexec -n 4 ./CNS3d.gnu.MPI.OMP.ex $inFile > $outFile
mpiexec -n 4 ./CNS3d.gnu.MPI.OMP.ex $inFile

