#!/bin/bash
#SBATCH -n 64
#SBATCH -p normal
#SBATCH -t 2:00:00
source options.sh
cwd=$(pwd)
fem="$cwd/../FEM/bspfem"
meshdir="$cwd/../Meshes"
for mesh in "${mesh_arr[@]}"
do
for poly in $meshdir/$mesh*.m
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
done
