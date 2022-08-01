gfortran -fPIC -fopenmp -w -O3 -std=legacy test_crn.f ../src/qerrfun.f
./a.out
python plot_curv.py
