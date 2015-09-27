![$\Omega_h$](omega_h.png?raw=true "Omega sub h")

# Compiling

Copy one of the *config.mk files and call
it config.mk, or write your own.
Then you can run `make`.

# Included programs:

warp.exe: Runs a 2D mesh adaptation test
including mesh motion with a size field based
on the hessian of a dye field.

warp_3d.exe: Runs a 3D mesh adaptation test
similar to test_warp.exe except the size field
is constant.

vtk_surfer.exe: A VTK file post-processing tool,
run with no arguments for more information.

vtkdiff.exe: A VTK file comparison tool with
floating point tolerances, run with no arguments
or with -help for information.
