#!/bin/bash

#mydir=./test32_time1
#mydir=./test64_time1
#mydir=./test128_time1
#mydir=./test256_time1

mydir=./test64_time1_AMR1

mkdir output

make

mpiexec -n 4 ./main2d.gnu.MPI.ex ./inputs

tree output

ls -1 plt*/Header | tee movie.visit

rm -rf $mydir
mkdir $mydir
cp inputs ./$mydir/inputs
mv movie.visit ./$mydir
mv plt* ./$mydir
mv output ./$mydir