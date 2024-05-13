#!/bin/env bash
# recompile the file if needed
make
maxVal=100

for (( n=10; n<=$maxVal; n+=10 ))
do
    echo "n=$n"
    echo "Running hqr 10 times"
    for (( l=1; l<=10; l+=1 ))
    do
        ./timeHSEQR.exe -n $n
    done
done
