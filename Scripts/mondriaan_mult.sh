#!/bin/bash
source options
for proc in "${parr[@]}"
do
  ./mondriaan.sh $1 $proc
done
