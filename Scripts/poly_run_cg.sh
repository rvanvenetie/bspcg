#!/bin/bash
#SBATCH -p short
#SBATCH -n 64
#SBATCH -t 30:00
source options.sh
cwd=$(pwd)
cg="$cwd/../Parallel/main"
meshdir="$cwd/../Meshes"
for poly in $meshdir/poly*.m
do
  for p in "${parr[@]}"
  do
    polyname="$(basename "$poly" .m)"
    matname="$meshdir/$polyname/$polyname"
		if [ -f "$matname.m-P$p" ] && [ -f "$matname.mtx-P$p" ]
		then
		  $runcmd $cg $matname.mtx $p
		fi
  done
done
