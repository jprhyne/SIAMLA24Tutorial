#!/bin/env bash
# recompile the file if needed
make
maxVal=1000

for (( n=100; n<=$maxVal; n+=100 ))
do
    echo "n=$n"
    echo "Running hqr 10 times"
    for (( l=1; l<=10; l+=1 ))
    do
        ./experiments.exe -n $n
    done
done
