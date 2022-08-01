gfortran -w -fopenmp -O3 -std=legacy test_quad.f ../src/get_geoms.f -L/mnt/home/mrachh/lib -lfmm2dbie -lfmm2d -lopenblas -lgomp
./a.out
#valgrind ./a.out
