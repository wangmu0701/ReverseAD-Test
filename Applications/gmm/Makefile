ReverseADHOME = $(HOME)/packages/reversead
ADOLCHOME = $(HOME)/packages/adolc_edge

CXX=g++

all : gmm_reversead gmm_adolc gmm_fullhess

gmm_reversead : gmm.cpp gmm.h defs.h utils.h
	$(CXX) -std=c++11 -O3 -D USING_REVERSEAD -I$(ReverseADHOME)/include gmm.cpp -o $@ -L$(ReverseADHOME)/lib -lreversead

gmm_adolc : gmm.cpp gmm.h defs.h utils.h
	$(CXX) -std=c++11 -O3 -D USING_ADOLC -I$(ADOLCHOME)/include gmm.cpp -o $@ -L$(ADOLCHOME)/lib64 -ladolc

gmm_fullhess : gmm.cpp gmm.h defs.h utils.h
	$(CXX) -std=c++11 -O3 -D USING_ADOLC -D ADOLC_FULLHESS -I$(ADOLCHOME)/include gmm.cpp -o $@ -L$(ADOLCHOME)/lib64 -ladolc

clean:
	rm gmm_reversead gmm_adolc gmm_fullhess

.PHONY:	clean
