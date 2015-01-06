#!/bin/bash
source options.sh
for proc in "${parr[@]}"
do
  ./mondriaan.sh $1 $proc
done
