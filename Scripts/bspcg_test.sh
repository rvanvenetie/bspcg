#!/bin/bash
#SBATCH -n 64
#SBATCH -t 30:00
source options.sh
cwd=$(pwd)
bspcg="$cwd/../Parallel/main"
matdir="$cwd/../mtxMatrices"
for n in "${narr[@]}"
do
  for sp in "${sparr[@]}"
	do
	  for p in "${parr[@]}"
		do
			matname="randmat_${n}_${sp/./_}"
			matfile="$matdir/$matname/$matname.mtx"
			$runcmd $bspcg $p $matfile $loadb $kmax $eps $time $debug 
		done
	done
done
