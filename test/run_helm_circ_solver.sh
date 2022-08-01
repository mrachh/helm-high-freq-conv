gfortran -fPIC -fopenmp -w -O3 -std=legacy test_helm_circ_solver.f ../src/get_geoms.f ../src/helm_galerkin_routs.f ../src/qerrfun.f -L/usr/local/lib -lfmm2dbie -lfmm2d -framework accelerate
./a.out
