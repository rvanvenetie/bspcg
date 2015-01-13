#!/bin/bash
source options.sh
cwd=$(pwd)
fem="$cwd/../FEM/bspfem"
meshdir="$cwd/../Meshes"
for poly in $meshdir/poly*.m
do
  for p in "${parr[@]}"
  do
    polyname= $(basename "$poly" .m)
    matname="$meshdir/$polyname/$polyname.m"
    $fem $matname $p
  done
done
