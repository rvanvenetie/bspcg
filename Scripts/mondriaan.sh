#!/bin/bash
cwd=$(pwd)
mond="$cwd/../Mondriaan4/tools/Mondriaan"
matdir="$cwd/../Meshes/$1"
matname="$matdir/$1.mtx"
mondopts="0.03 -SquareMatrix_DistributeVectorsEqual=yes"
echo "Generating mondriaan file for $1 with $2 processors"
$mond $matname $2 $mondopts
