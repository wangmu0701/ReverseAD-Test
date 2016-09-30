CXX=g++
CXXFLAGS=-std=c++11

REVERSEADPATH=$(HOME)/packages/reversead/

all : stream_vel stream_vel_reversead gradient

stream_vel : stream_vel.cpp stream_vel.hpp
	$(CXX) $(CXXFLAGS) -o $@ stream_vel.cpp

stream_vel_reversead : stream_vel.cpp stream_vel.hpp
	$(CXX) $(CXXFLAGS) -D USING_REVERSEAD -o $@ stream_vel.cpp -I$(REVERSEADPATH)/include -L$(REVERSEADPATH)/lib -lreversead

gradient : stream_vel.cpp stream_vel.hpp
	$(CXX) $(CXXFLAGS) -D USING_REVERSEAD -D CHECK_GRADIENT -o $@ stream_vel.cpp -I$(REVERSEADPATH)/include -L$(REVERSEADPATH)/lib -lreversead

clean:
	rm -f stream_vel stream_vel_reversead gradient

.PHONY: clean
