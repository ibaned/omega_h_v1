#!/bin/bash -ex

$VALGRIND ./bin/loop.exe
if [ "$LOOP_MODE" = "cuda" ]; then
  return
fi
$VALGRIND ./bin/box.exe --file ./data/box.vtu --dim 2 --refinements 6
$VALGRIND ./bin/vtkdiff.exe ./data/box.vtu ./data/box.vtu
$VALGRIND ./bin/node_ele.exe ./data/xgc.node ./data/xgc.ele ./data/xgc.vtu
$VALGRIND ./bin/from_gmsh.exe ./data/cube.msh ./data/cube.vtu
if [ "$USE_MPI" = "1" ]; then
  $MPIRUN -np 2 $VALGRIND ./bin/migrate.exe
  $MPIRUN -np 2 $VALGRIND ./bin/conform.exe
  $MPIRUN -np 2 $VALGRIND ./bin/partition.exe ./data/box.vtu ./data/split.pvtu
  $MPIRUN -np 2 $VALGRIND ./bin/one_refine.exe ./data/split.pvtu ./data/one_ref.pvtu
  $MPIRUN -np 2 $VALGRIND ./bin/one_refine.exe ./data/one_ref.pvtu ./data/two_ref.pvtu
  $MPIRUN -np 2 $VALGRIND ./bin/one_coarsen.exe ./data/split.pvtu ./data/one_cor.pvtu
  $MPIRUN -np 2 $VALGRIND ./bin/one_coarsen.exe ./data/one_cor.pvtu ./data/two_cor.pvtu
fi
$VALGRIND ./bin/warp.exe
if [ "$PATIENT" = "1" ]; then
  $VALGRIND ./bin/warp_3d.exe
fi
