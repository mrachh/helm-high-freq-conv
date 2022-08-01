gfortran -w -O3 -fopenmp -std=legacy test_quad.f ../src/get_geoms.f ../src/helm_galerkin_routs.f -L/usr/local/lib -lfmm2dbie -lfmm2d -framework accelerate
./a.out
