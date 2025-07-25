rm ./a.out
gfortran -fopenmp -fPIC -w -O3 -std=legacy diamond_example.f ../src/get_geoms.f ../src/qerrfun.f ../src/helm_galerkin_routs_paper.f -L/usr/local/lib -lfmm2dbie -lfmm2d -lopenblas -lgomp
./a.out | tee results/results_jul23_2025.dat
