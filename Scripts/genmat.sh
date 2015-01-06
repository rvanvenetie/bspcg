#!/bin/bash
cwd=$(pwd)
genmat="$(pwd)/../Genmat/rand_sym_mat"
matdir="../mtxMatrices/"
matname="randmat_$1_${2/./_}"
cd $matdir
mkdir -p $matname
cd $matname
echo "Generating random matrix with n=$1 and sparsity=$2"
$genmat $1 $2 > "$matname.mtx"
cd $cwd
./mondriaan_mult.sh $matname
