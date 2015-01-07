#!/bin/bash
# Use with matrixname as parameter
source options.sh
for proc in "${parr[@]}"
do
  ./mondriaan.sh $1 $proc
done
