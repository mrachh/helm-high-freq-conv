gfortran -w -O3 -std=legacy test_quad.f -L/mnt/home/mrachh/lib -lfmm2dbie -lfmm2d -lopenblas
./a.out
#valgrind ./a.out
