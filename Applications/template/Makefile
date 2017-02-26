INSTANCE_PATH=./mesh
#INSTANCE_PATH=./gmm
#INSTANCE_PATH=./hand

REVERSEADPATH=$(HOME)/packages/reversead
ADOLCPATH=$(HOME)/packages/adolc
LIBTAYLORPATH=$(HOME)/packages/libtaylor

DUMMYSCALAR=-D DUMMY_SCALAR

DEVORDER=2
ADOLCMETHOD=3
REVERSEADMETHOD=2
LIBTAYLORMETHOD=0

ADOLCFLAGS=-D ORDER=$(DEVORDER) -D METHOD=$(ADOLCMETHOD)
REVERSEADFLAGS=-D ORDER=$(DEVORDER) -D METHOD=$(REVERSEADMETHOD)
LIBTAYLORFLAGS=-D ORDER=$(DEVORDER) -D METHOD=$(LIBTAYLORMETHOD)

CXX=g++

CXXFLAGS=-std=c++11 -O3 -I$(INSTANCE_PATH)

all : test_reversead test_reversead_d test_adolc test_adolc_d test_libtaylor

test_reversead : func_eval.o test_reversead.o
	$(CXX) $(CXXFLAGS) $(REVERSEADFLAGS) -o $@ $^ -I$(INSTANCE_PATH) -I$(REVERSEADPATH)/include -L$(REVERSEADPATH)/lib -lreversead

test_reversead.o : test_reversead.cpp
	$(CXX) $(CXXFLAGS) $(REVERSEADFLAGS) -o $@ -c $^ -I$(REVERSEADPATH)/include

test_reversead_d : func_eval.o test_reversead_d.o
	$(CXX) $(CXXFLAGS) $(REVERSEADFLAGS) $(DUMMYSCALAR) -o $@ $^ -I$(INSTANCE_PATH) -I$(REVERSEADPATH)/include -L$(REVERSEADPATH)/lib -lreversead

test_reversead_d.o : test_reversead.cpp
	$(CXX) $(CXXFLAGS) $(REVERSEADFLAGS) $(DUMMYSCALAR) -o $@ -c $^ -I$(REVERSEADPATH)/include

test_adolc : func_eval.o test_adolc.o
	$(CXX) $(CXXFLAGS) $(ADOLCFLAGS) -o $@ $^ -I$(INSTANCE_PATH) -I$(ADOLCPATH)/include -L$(ADOLCPATH)/lib64 -ladolc

test_adolc.o : test_adolc.cpp
	$(CXX) $(CXXFLAGS) $(ADOLCFLAGS) -o $@ -c $^ -I$(ADOLCPATH)/include

test_adolc_d : func_eval.o test_adolc_d.o
	$(CXX) $(CXXFLAGS) $(ADOLCFLAGS) $(DUMMYSCALAR) -o $@ $^ -I$(INSTANCE_PATH) -I$(ADOLCPATH)/include -L$(ADOLCPATH)/lib64 -ladolc

test_adolc_d.o : test_adolc.cpp
	$(CXX) $(CXXFLAGS) $(ADOLCFLAGS) $(DUMMYSCALAR) -o $@ -c $^ -I$(ADOLCPATH)/include

func_eval.o : $(INSTANCE_PATH)/func_eval.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $^ -I$(INSTANCE_PATH)


test_libtaylor : func_eval.o test_libtaylor.o
	$(CXX) $(CXXFLAGS) $(LIBTAYLORFLAGS) -o $@ $^ -I$(INSTANCE_PATH) -I$(LIBTAYLORPATH)/include

test_libtaylor.o : test_libtaylor.cpp
	$(CXX) $(CXXFLAGS) $(LIBTAYLORFLAGS) -o $@ -c $^ -I$(LIBTAYLORPATH)/include



.PHONY : clean

clean :
	rm test_reversead test_adolc test_reversead_d test_adolc_d test_libtaylor *.o
