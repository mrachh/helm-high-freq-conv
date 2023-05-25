rm ./a.out
gfortran -fopenmp -fPIC -w -O3 -std=legacy test_helm_cavity_solver_new_proj_lin.f ../src/get_geoms.f ../src/qerrfun.f ../src/helm_galerkin_routs.f -L/mnt/home/mrachh/lib -lfmm2dbie -lfmm2d -lopenblas -lgomp
./a.out | tee results_sep30.dat
