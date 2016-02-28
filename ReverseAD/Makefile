ReverseadHome = $(HOME)/packages/reversead
CXX=/usr/local/bin/g++

all : reversead

reversead : reversead.cpp ./../synthetic.hpp
	$(CXX) -std=c++11 -I$(ReverseadHome)/include reversead.cpp -o $@ -L$(ReverseadHome)/lib -lreversead

clean :
	rm reversead
