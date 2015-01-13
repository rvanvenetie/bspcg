#!/bin/bash
source options.sh
./m2fem.sh $1 
for proc in "${parr[@]}"
do
	./m2dm.sh $1 $proc
	./mondriaan.sh $1 $proc
done
