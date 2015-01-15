#!/bin/bash
mesh_arr=("poly5")
#mesh_arr=("lshaped_" "poly5" "poly8")
if [[ `whoami` == "bissstud" ]]; then
	cd $HOME/Students14/Jan/bspcg/Scripts/
	echo "Running on Cartesius $@"
	runcmd="srun"
else
	echo "Running locally $@"
	runcmd=""
fi
echo -e "\n\033[1;31mtime\t`date`\033[0;39m"

narr=(16 32 64)
sparr=(0.01)
parr=(1 2 4)
narr=(64 128 256 512 1024 2048 4096)
sparr=(0.01 0.02 0.03 0.04 0.05 0.1 0.2)
parr=(1 2 4 8 16 32 64)
kmax=100000
eps=0.01
time=1
debug=0
loadb=0
