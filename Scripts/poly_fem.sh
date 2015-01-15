#!/bin/bash
for poly in ../Meshes/$1*.m
do
  ./fem_mult.sh $(basename "$poly" .m)
done
