#!/bin/bash
for poly in ../Meshes/poly*.m
do
  ./fem_mult.sh $(basename "$poly" .m)
done
