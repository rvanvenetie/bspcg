#!/bin/bash
procs=(1 2 3 4 5);
for proc in "${procs[@]}"
do
  ./mondriaan.sh $1 $proc
done
