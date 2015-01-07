#!/bin/bash
#SBATCH -n 64
#SBATCH -t 30:00
runcmd="srun"
cwd=/home/bissstud/Students14/Jan/bspcg/Scripts
bspcg="$cwd/../Parallel/main"
matdir="$cwd/../mtxMatrices"
source options.sh
for n in "${narr[@]}"
do
  for sp in "${sparr[@]}"
	do
	  for p in "${parr[@]}"
		do
			matname="randmat_${n}_${sp/./_}"
			matfile="$matdir/$matname/$matname.mtx"
			$runcmd $bspcg $matfile $loadb $p $kmax $eps $time $debug 
		done
	done
done
