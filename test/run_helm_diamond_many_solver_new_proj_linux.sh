rm ./a.out
gfortran -fopenmp -fPIC -w -O3 -std=legacy test_helm_diamond_many_solver_new_proj.f ../src/get_geoms.f ../src/qerrfun.f ../src/helm_galerkin_routs.f -L/mnt/home/mrachh/lib -lfmm2dbie -lfmm2d -lopenblas -lgomp
./a.out | tee results_sep29.dat
