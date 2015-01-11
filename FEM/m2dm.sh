cwd=$(pwd)
mond="$cwd/../Mondriaan4/tools/Mondriaan"
m2mtx="$cwd/../FEM/m2mtx"
mtx2dm="$cwd/../FEM/mtx2dm"
mesh=$1
P=$2
$m2mtx $mesh.m > $mesh.mtx
$mond $mesh.mtx $2 0.1
$mtx2dm $mesh.m $2 > $mesh.m-P$2
