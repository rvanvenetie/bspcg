#!/bin/bash
narr=(16 32 64 128 256 512 1024 2048 4096)
sparr=(0.01 0.02 0.03 0.04 0.05 0.1 0.2)
for n in "${narr[@]}"
do
  for sp in "${sparr[@]}"
	do
    ./genmat.sh $n $sp
	done
done
