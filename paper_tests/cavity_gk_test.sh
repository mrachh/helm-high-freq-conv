rm ./a.out
gfortran -fopenmp -fPIC -w -O3 -std=legacy cavity_gk_analytic.f ../src/get_geoms.f ../src/qerrfun.f ../src/helm_galerkin_routs_paper.f -L/usr/local/lib -lfmm2dbie -lfmm2d -lopenblas -lgomp
./a.out | tee results_sep30.dat
