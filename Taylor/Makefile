all: taylorTest

taylorTest: taylorTest.cpp ./../synthetic.hpp
	/usr/local/bin/g++ -g -I/Users/muwang/packages/adolc_edge/include/ -I/usr/local/include/ -o taylorTest taylorTest.cpp -L/Users/muwang/packages/adolc_edge/lib64 -ladolc

clean:
	rm taylorTest
