gfortran -fPIC -w -O3 -std=legacy test_helm_circ_solver_paper.f ../src/get_geoms.f ../src/qerrfun.f ../src/helm_galerkin_routs.f -L/mnt/home/mrachh/lib -lfmm2dbie -lfmm2d -lopenblas
./a.out
