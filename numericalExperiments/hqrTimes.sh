#!/bin/env bash
# recompile if needed
make
maxVal=100

for (( n=10; n<=$maxVal; n+=10 ))
do
    echo "n=$n"
    echo "Running fortran hqr 10 times"
    for (( l=1; l<=10; l+=1 ))
    do
        echo -n "TQE:"
        ./hqr/test_hqr2schur_fortran.exe -n $n -e
    done
done
maxVal=1000

for (( n=100; n<=$maxVal; n+=100 ))
do
    echo "n=$n"
    echo "Running fortran hqr 10 times"
    for (( l=1; l<=10; l+=1 ))
    do
        echo -n "TQE:"
        ./hqr/test_hqr2schur_fortran.exe -n $n -e
    done
done
