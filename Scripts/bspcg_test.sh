#!/bin/bash
cwd=$(pwd)
bspcg="$cwd/../Parallel/main"
matdir="$cwd/../mtxMatrices"
source options.sh
for n in "${narr[@]}"
do
  for sp in "${sparr[@]}"
	do
	  ./genmat.sh $n $sp
	  for p in "${parr[@]}"
		do
			matname="randmat_${n}_${sp/./_}"
			matfile="$matdir/$matname/$matname.mtx"
			$bspcg $matfile $loadb $p $kmax $eps $time $debug
		done
	done
done
