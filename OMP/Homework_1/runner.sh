#!/bin/bash
for p in 1 2 4 8 16; do
    for ((N=1000; N <= 30000; N += 1000)) do
        ./a.out 0 100 50 "$N" 0.5 "$p" 
    done
    # ./a.out 0 100 50 1000000 0.5 "$p"
done