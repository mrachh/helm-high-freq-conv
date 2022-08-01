gfortran -fPIC -w -O3 -std=legacy test_helm_circ_solver.f ../src/get_geoms.f ../src/helm_galerkin_routs.f -L/usr/local/lib -lfmm2dbie -lfmm2d -framework accelerate
./a.out
