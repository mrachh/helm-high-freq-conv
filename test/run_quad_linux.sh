gfortran -w -fopenmp -O3 -std=legacy test_quad.f ../src/get_geoms.f ../src/helm_galerkin_routs.f ../src/qerrfun.f -L/mnt/home/mrachh/lib -lfmm2dbie -lfmm2d -lblas -llapack -lgomp
./a.out
#valgrind ./a.out
