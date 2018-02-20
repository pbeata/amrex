#!/bin/bash


# Base input file
prefix='inputs'
fileType='amr'
outType='out'
outDir='output'
exe='CNS3d.gnu.MPI.ex'

# AMReX Distribution Mapping Strategies:
# SF = Space Filling Curve
# KS = Knapsack
# RR = Round Robin
DistMapStrats='SF KS RR'

# Simulation time in the test:
# full = 1.2e-5 seconds
# half = 0.6e-6 seconds
SimTime='half full'
SimTime='half'

# MPI processes:
MPI='1 2 4 8 16 32 40 80'
MPI='8'



# LOOP OVER THE TEST CONDITIONS
for strat in $DistMapStrats
do
  file1=$prefix'_'$strat
  for time in $SimTime
  do
    file2='_'$time'.'$fileType
    inFile=$file1$file2
    echo -e '\n'$inFile
    for N in $MPI
    do
      outTest=$outDir'_'$strat'_'$time'_'$N'MPI'
      outFile=$strat'_'$time'_'$N'MPI.'$outType
      
      # runTest='mpiexec -n '$N' ./CNS3d.gnu.MPI.ex ./'$inFile' > '$outFile
      # echo $runTest
      # echo $outTest
      
      #==========================================
      # RUN THE TEST
      mkdir output
      # mpiexec -n $N ./$exe ./$inFile > $outFile
      mv output $outTest
      mv $outFile $outTest/$outFile 
      #==========================================

    done
  done
done

