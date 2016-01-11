#!/bin/bash -ex

$VALGRIND ./bin/print_swap_edges.exe

$VALGRIND ./bin/box.exe --file ./data/box.vtu --dim 2 --refinements 6

$VALGRIND ./bin/node_ele.exe ./data/xgc.node ./data/xgc.ele ./data/xgc.vtu

#$VALGRIND ./bin/warp.exe

#$VALGRIND ./bin/warp_3d.exe

$MPIRUN -np 2 $VALGRIND ./bin/migrate.exe

$MPIRUN -np 2 $VALGRIND ./bin/conform.exe
