#!/bin/bash
cwd=$(pwd)
m2fem="$cwd/../MeshTools/m2fem"
meshdir="$cwd/../Meshes/"
meshname=$1
# Go to the mesh directory.
# Make a directory with the name of the current mesh file
# Mesh should be in the root of the mesh dir
# Copy the mesh file to the newly created dir
# Create the fem matrix in this directy
cd $meshdir
mkdir -p $meshname
#Copy file into newly created dir
cp $meshname.m $meshname/$meshname.m
#Got to newly created dir
cd $meshname
echo "Generating fem system for $meshname"
#Generate hypegraph
$m2fem $meshname.m > $meshname.mtx 2> $meshname.mtx-b
