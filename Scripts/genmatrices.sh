#!/bin/bash
source options.sh
for n in "${narr[@]}"
do
  for sp in "${sparr[@]}"
	do
    ./genmat.sh $n $sp
	done
done
