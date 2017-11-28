#!/bin/bash
for P in 1 2 4 8 16; do
    for ((N=10000; N <= 30000000; N += 1000000)) do
        let "M = (N + P - 1) / P"
        ./a.out "$N" "$M" "$P"
        #for ((M=10000; M <= N; M += 1000000)) do
        #    ./a.out "$N" "$M" "$P"
        #done
    done
    echo "$P done"
done
