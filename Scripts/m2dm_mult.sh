#!/bin/bash
source options.sh
for proc in "${parr[@]}$"
do
	./m2dm.sh $1 $proc
done
