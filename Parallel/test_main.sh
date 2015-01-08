#!/bin/bash
cwd=$(pwd)
bspcg="$cwd/main"
matdir="$cwd/../mtxMatrices"
matfile="$matdir/$1/$1.mtx"

$bspcg $matfile
