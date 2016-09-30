CXX=g++
CXXFLAGS=-std=c++11
all : stream_vel

stream_vel : driver.o stream_vel_timedep.o
	$(CXX) $(CXXFLAGS) -o $@ $^ 


%.o : %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -f *.o stream_vel

.PHONY: clean
