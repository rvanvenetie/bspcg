#!/bin/bash
cwd=$(pwd)
mond="$cwd/../Mondriaan4/tools/Mondriaan"
m2mtx="$cwd/../MeshTools/m2mtx"
mtx2dm="$cwd/../MeshTools/mtx2dm"
meshdir="$cwd/../Meshes/"
meshname=$1
p=$2
# Go to the mesh directory.
# Make a directory with the name of the current mesh file
# Mesh should be in the root of the mesh dir
# Copy the mesh file to the newly created dir
# Create a distributed mesh in this directory for proc given in $2
cd $meshdir
mkdir -p $meshname
#Copy file into newly created dir
cp $meshname.m $meshname/$meshname.m
#Got to newly created dir
cd $meshname
echo "Generating distributed mesh for $meshname with p=$p"
#Generate hypegraph
$m2mtx $meshname.m > $meshname.mtx_tmp
#Generate distributed hypergraph
$mond $meshname.mtx_tmp $p 0.03
#Generate distributed mesh
$mtx2dm $meshname.m $p > $meshname.m-P$p
rm *.mtx_tmp*
