#!/bin/bash
source options.sh
for proc in "${parr[@]}$"
do
	./m2dm.sh $1 $proc
	./m2fem.sh $1 $proc
	./mondriaan_mult.sh $1 $proc
done
