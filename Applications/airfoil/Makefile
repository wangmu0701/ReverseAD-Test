ReverseADHOME = $(HOME)/packages/reversead
ADOLCHOME = $(HOME)/packages/adolc_edge

CXX=/usr/local/bin/g++

all : airfoil airfoil_revad airfoil_deriv airfoil_adolc airfoil_adolc_dev

airfoil : airfoil_r.cpp input_r.cpp routines_r.hpp const.h
	$(CXX) -std=c++11 -O3 -I$(ReverseADHOME)/include airfoil_r.cpp input_r.cpp -o $@ -L$(ReverseADHOME)/lib -lreversead

airfoil_revad : airfoil_r.cpp input_r.cpp routines_r.hpp const.h
	$(CXX) -std=c++11 -O3 -D USING_REVERSEAD -I$(ReverseADHOME)/include airfoil_r.cpp input_r.cpp -o $@ -L$(ReverseADHOME)/lib -lreversead

airfoil_deriv : airfoil_r.cpp input_r.cpp routines_r.hpp const.h
	$(CXX) -std=c++11 -O3 -D USING_REVERSEAD -D COMPUTE_DERIVATIVE -I$(ReverseADHOME)/include airfoil_r.cpp input_r.cpp -o $@ -L$(ReverseADHOME)/lib -lreversead

airfoil_adolc : airfoil_adolc.cpp input_r.cpp routines_r.hpp const.h
	$(CXX) -std=c++11 -O3 -I$(ADOLCHOME)/include airfoil_adolc.cpp input_r.cpp -o $@ -L$(ADOLCHOME)/lib64 -ladolc

airfoil_adolc_dev : airfoil_adolc.cpp input_r.cpp routines_r.hpp const.h
	$(CXX) -std=c++11 -O3 -D COMPUTE_DERIVATIVE -I$(ADOLCHOME)/include airfoil_adolc.cpp input_r.cpp -o $@ -L$(ADOLCHOME)/lib64 -ladolc

clean :
	rm airfoil airfoil_revad airfoil_deriv airfoil_adolc airfoil_adolc_dev
