#!/bin/bash
source options.sh
cwd=$(pwd)
cg="$cwd/../Parallel/main"
meshdir="$cwd/../Meshes"
for poly in $meshdir/poly*.m
do
  for p in "${parr[@]}"
  do
    polyname="$(basename "$poly" .m)"
    matname="$meshdir/$polyname/$polyname.m"
    $cg $matname $p
  done
done
