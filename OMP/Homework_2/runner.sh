#!/bin/bash
for P in 1 2 4 8 16; do
    for ((N=10; N <= 30; N += 10)) do
        ./a.out "$N" "$N" "$P"
        # for ((M=10000; M <= N; M += 10000)) do
        #     ./a.out "$N" "$M" "$P"
        # done
    done
    # echo "$P done"
done
