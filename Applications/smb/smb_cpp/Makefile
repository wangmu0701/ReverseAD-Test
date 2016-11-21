CXX=g++
CXXFLAGS=-O3 --std=c++11
REVERSEADPATH=$(HOME)/packages/reversead
ADOLCPATH=$(HOME)/packages/adolc
LIBTAYLORPATH=$(HOME)/packages/libtaylor

ADTOOL=-D USING_ADTOOL
USEREVERSEAD=$(ADTOOL) -D USING_REVERSEAD
USEADOLC=$(ADTOOL) -D USING_ADOLC
REV1=-D COMPUTE_GRADIENT #-D CHECK_GRADIENT
REV2=-D COMPUTE_HESSIAN #-D CHECK_HESSIAN
REV3=-D COMPUTE_THIRD #-D CHECK_THIRD
REVH=-D COMPUTE_HIGHER_ORDER #-D CHECK_HIGHER_ORDER
#TAYLORSPARSE=-D USING_SPARSITY

all : smb_main rev_gradient rev_hessian rev_third rev_high\
      adolc_gradient adolc_hessian adolc_high\
      taylor_gradient taylor_hessian taylor_third
#all : adolc_hessian

smb_main : smb_main.cpp eval_func_chem.hpp
	$(CXX) $(CXXFLAGS) -o $@ smb_main.cpp

rev_gradient : smb_main.cpp eval_func_chem.hpp
	$(CXX) $(CXXFLAGS) $(USEREVERSEAD) $(REV1) -o $@ smb_main.cpp -I$(REVERSEADPATH)/include -L$(REVERSEADPATH)/lib -lreversead

rev_hessian : smb_main.cpp eval_func_chem.hpp
	$(CXX) $(CXXFLAGS) $(USEREVERSEAD) $(REV2) -o $@ smb_main.cpp -I$(REVERSEADPATH)/include -L$(REVERSEADPATH)/lib -lreversead

rev_third : smb_main.cpp eval_func_chem.hpp
	$(CXX) $(CXXFLAGS) $(USEREVERSEAD) $(REV3) -o $@ smb_main.cpp -I$(REVERSEADPATH)/include -L$(REVERSEADPATH)/lib -lreversead

rev_high : smb_main.cpp eval_func_chem.hpp
	$(CXX) $(CXXFLAGS) $(USEREVERSEAD) $(REVH) -o $@ smb_main.cpp -I$(REVERSEADPATH)/include -L$(REVERSEADPATH)/lib -lreversead

adolc_gradient : smb_main.cpp eval_func_chem.hpp
	$(CXX) $(CXXFLAGS) $(USEADOLC) $(REV1) -o $@ smb_main.cpp -I$(ADOLCPATH)/include -L$(ADOLCPATH)/lib64 -ladolc

adolc_hessian : smb_main.cpp eval_func_chem.hpp
	$(CXX) $(CXXFLAGS) $(USEADOLC) $(REV2) -o $@ smb_main.cpp -I$(ADOLCPATH)/include -L$(ADOLCPATH)/lib64 -ladolc

adolc_high : smb_main.cpp eval_func_chem.hpp
	$(CXX) $(CXXFLAGS) $(USEADOLC) $(REVH) -o $@ smb_main.cpp -I$(ADOLCPATH)/include -L$(ADOLCPATH)/lib64 -ladolc

taylor_gradient : smb_libtaylor.cpp eval_func_chem.hpp
	$(CXX) $(CXXFLAGS) $(REV1) -o $@ smb_libtaylor.cpp -I$(LIBTAYLORPATH)/include

taylor_hessian : smb_libtaylor.cpp eval_func_chem.hpp
	$(CXX) $(CXXFLAGS) $(REV2) -o $@ smb_libtaylor.cpp -I$(LIBTAYLORPATH)/include

taylor_third : smb_libtaylor.cpp eval_func_chem.hpp
	$(CXX) $(CXXFLAGS) $(REV3) -o $@ smb_libtaylor.cpp -I$(LIBTAYLORPATH)/include

clean:
	rm smb_main rev_gradient rev_hessian rev_third rev_high\
           adolc_gradient adolc_hessian adolc_high\
           taylor_gradient taylor_hessian taylor_third

.PHONY : clean
