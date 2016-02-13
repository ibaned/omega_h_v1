#!/bin/bash -ex

$VALGRIND ./bin/loop.exe
if [ "$LOOP_MODE" = "cuda" ]; then
  return
fi
$VALGRIND ./bin/box.exe --file ./scratch/box.vtu --dim 2 --refinements 6
$VALGRIND ./bin/vtkdiff.exe ./scratch/box.vtu ./scratch/box.vtu
$VALGRIND ./bin/node_ele.exe ./data/xgc.node ./data/xgc.ele ./scratch/xgc.vtu
$VALGRIND ./bin/from_gmsh.exe ./data/cube.msh ./scratch/cube.vtu
if [ "$USE_MPI" = "1" ]; then
  $MPIRUN -np 2 $VALGRIND ./bin/migrate.exe ./scratch
  $MPIRUN -np 2 $VALGRIND ./bin/conform.exe ./scratch
  $MPIRUN -np 2 $VALGRIND ./bin/partition.exe ./scratch/box.vtu ./scratch/split.pvtu
  $MPIRUN -np 2 $VALGRIND ./bin/one_refine.exe ./scratch/split.pvtu ./scratch/one_ref.pvtu
  $MPIRUN -np 2 $VALGRIND ./bin/one_refine.exe ./scratch/one_ref.pvtu ./scratch/two_ref.pvtu
  $MPIRUN -np 2 $VALGRIND ./bin/one_coarsen.exe ./scratch/split.pvtu ./scratch/one_cor.pvtu
  $MPIRUN -np 2 $VALGRIND ./bin/one_coarsen.exe ./scratch/one_cor.pvtu ./scratch/two_cor.pvtu
fi
$VALGRIND ./bin/warp.exe ./scratch
if [ -e ./gold/2d_warp_0008.vtu ]; then
  $VALGRIND ./bin/vtkdiff.exe ./scratch/warp_0008.vtu ./gold/2d_warp_0008.vtu
  diff ./scratch/warp_0008.vtu ./gold/2d_warp_0008.vtu
fi
cp ./scratch/warp_0008.vtu ./gold/2d_warp_0008.vtu
if [ "$PATIENT" = "1" ]; then
  $VALGRIND ./bin/warp_3d.exe ./scratch
  if [ -e ./gold/warp_0016.vtu ]; then
    $VALGRIND ./bin/vtkdiff.exe ./scratch/warp_0016.vtu ./gold/warp_0016.vtu
    diff ./scratch/warp_0016.vtu ./gold/warp_0016.vtu
  fi
  cp ./scratch/warp_0016.vtu ./gold/warp_0016.vtu
fi
