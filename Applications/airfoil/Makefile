ReverseADHOME = $(HOME)/packages/reversead
ADOLCHOME = $(HOME)/packages/adolc

#CXX=/usr/local/bin/g++
CXX=g++

all : airfoil airfoil_revad airfoil_deriv airfoil_adolc airfoil_adolc_dev check_hess

airfoil : airfoil_r.cpp airfoil_func.hpp input_r.cpp routines_r.hpp const.h
	$(CXX) -std=c++11 -g -O0 -I$(ReverseADHOME)/include airfoil_r.cpp input_r.cpp -o $@ -L$(ReverseADHOME)/lib -lreversead

airfoil_revad : airfoil_r.cpp airfoil_func.hpp input_r.cpp routines_r.hpp const.h
	$(CXX) -std=c++11 -O3 -D USING_REVERSEAD -I$(ReverseADHOME)/include airfoil_r.cpp input_r.cpp -o $@ -L$(ReverseADHOME)/lib -lreversead

airfoil_deriv : airfoil_r.cpp airfoil_func.hpp input_r.cpp routines_r.hpp const.h
	$(CXX) -std=c++11 -O3 -D USING_REVERSEAD -D COMPUTE_DERIVATIVE -I$(ReverseADHOME)/include airfoil_r.cpp input_r.cpp -o $@ -L$(ReverseADHOME)/lib -lreversead

airfoil_adolc : airfoil_adolc.cpp input_r.cpp routines_r.hpp const.h
	$(CXX) -std=c++11 -O3 -I$(ADOLCHOME)/include airfoil_adolc.cpp input_r.cpp -o $@ -L$(ADOLCHOME)/lib64 -ladolc

airfoil_adolc_dev : airfoil_adolc.cpp input_r.cpp routines_r.hpp const.h
	$(CXX) -std=c++11 -O3 -D COMPUTE_DERIVATIVE -I$(ADOLCHOME)/include airfoil_adolc.cpp input_r.cpp -o $@ -L$(ADOLCHOME)/lib64 -ladolc


check_hess : airfoil_r.cpp airfoil_func.hpp input_r.cpp routines_r.hpp const.h
	$(CXX) -std=c++11 -D HESSIAN_CHECK -I$(ReverseADHOME)/include airfoil_r.cpp input_r.cpp -o $@ -L$(ReverseADHOME)/lib -lreversead

clean :
	rm airfoil airfoil_revad airfoil_deriv airfoil_adolc airfoil_adolc_dev
