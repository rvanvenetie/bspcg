#!/bin/bash
#SBATCH -n 64
#SBATCH -p short
#SBATCH -t 30:00
source options.sh
cwd=$(pwd)
fem="$cwd/../FEM/bspfem"
meshdir="$cwd/../Meshes"
for poly in $meshdir/poly*.m
do
  for p in "${parr[@]}"
  do
    polyname="$(basename "$poly" .m)"
    matname="$meshdir/$polyname/$polyname"
		if [ -f "$matname.m-P$p" ] && [ -f "$matname.mtx-P$p" ]
		then
			$runcmd $fem $matname.m $p
		fi
  done
done
