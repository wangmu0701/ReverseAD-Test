all:    meshOpt hess_plain hess_reversead hess_adolc

REVERSEADPATH=$(HOME)/packages/reversead
ADOLCPATH=$(HOME)/packages/adolc

SOURCEFILES=./src/fcn3e_am.cpp ./src/tet.cpp ./src/mesh3.cpp ./src/Mesh.cpp ./src/opt3.cpp ./src/elem3.cpp ./src/pre3.cpp

meshOpt:
	g++ --std=c++11 -O3 -I./src -o $@ ./src/meshOpt.cpp $(SOURCEFILES)

hess_plain:
	g++ --std=c++11 -O3 -I./src -o $@ ./src/meshHess.cpp $(SOURCEFILES)
                        
hess_reversead:
	g++ --std=c++11 -O3 -D USING_REVERSEAD -I$(REVERSEADPATH)/include -I./src -o $@ ./src/meshHess.cpp $(SOURCEFILES) -L$(REVERSEADPATH)/lib -lreversead

hess_adolc:
	g++ --std=c++11 -O3 -D USING_ADOLC -I$(ADOLCPATH)/include -I./src -o $@ ./src/meshHess.cpp $(SOURCEFILES) -L$(ADOLCPATH)/lib64 -ladolc
	
.PHONY: meshOpt hess_plain hess_reversead hess_adolc

clean:
	rm *.tap meshOpt hess_* out.mesh
