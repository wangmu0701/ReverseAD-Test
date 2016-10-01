CXX=g++
CXXFLAGS=-std=c++11

REVERSEADPATH=$(HOME)/packages/reversead
ADOLCPATH=$(HOME)/packages/adolc

all : stream_vel stream_vel_reversead gradient hessian adolc

stream_vel : stream_vel.cpp stream_vel.hpp
	$(CXX) $(CXXFLAGS) -o $@ stream_vel.cpp

stream_vel_reversead : stream_vel.cpp stream_vel.hpp
	$(CXX) $(CXXFLAGS) -D USING_REVERSEAD -o $@ stream_vel.cpp -I$(REVERSEADPATH)/include -L$(REVERSEADPATH)/lib -lreversead

gradient : stream_vel.cpp stream_vel.hpp
	$(CXX) $(CXXFLAGS) -D USING_REVERSEAD -D REVERSEAD_GRADIENT -D CHECK_GRADIENT -o $@ stream_vel.cpp -I$(REVERSEADPATH)/include -L$(REVERSEADPATH)/lib -lreversead

hessian : stream_vel.cpp stream_vel.hpp
	$(CXX) $(CXXFLAGS) -D USING_REVERSEAD -D REVERSEAD_HESSIAN -D CHECK_HESSIAN -o $@ stream_vel.cpp -I$(REVERSEADPATH)/include -L$(REVERSEADPATH)/lib -lreversead

adolc : stream_vel.hpp adolc.cpp
	$(CXX) $(CXXFLAGS) -D USING_ADOLC -D CHECK_HESSIAN -o $@ adolc.cpp -I$(ADOLCPATH)/include -L$(ADOLCPATH)/lib64 -ladolc

clean:
	rm -f stream_vel stream_vel_reversead gradient hessian adolc

.PHONY: clean
