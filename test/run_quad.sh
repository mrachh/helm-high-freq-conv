gfortran -w -O3 -std=legacy test_quad.f ../src/get_geoms.f -L/usr/local/lib -lfmm2dbie -lfmm2d -framework accelerate
./a.out
