##########################################################
# This file is part of Rapsodia released under the LGPL. #
# The full COPYRIGHT notice can be found in the top      #
# level directory of the Rapsodia distribution           #
##########################################################
RAPSODIAHOME = $(HOME)/packages/Rapsodia
include $(RAPSODIAHOME)/MakeDefs.mk
include $(RAPSODIAHOME)/MakeRules.mk

default: driver
	./$^

GEN_DIR=RALib

RA_EXTRAS=$(RAPSODIAHOME)/hotCpp
IPATH+=-I$(GEN_DIR) -I$(RA_EXTRAS)

OBJS= \
$(addprefix $(RA_EXTRAS)/, $(addsuffix .o, $(HOTCPPNAMES))) \
driver.o

driver:  $(OBJS) $(GEN_DIR)/libRapsodia.a ../function.h
	$(CXX) $(CXXFLAGS) $(IPATH) -o $@ $^

$(GEN_DIR)/libRapsodia.a: FORCE
	$(RAPSODIAHOME)/Generator/generate.py -d 20 -s 20 -o 3 -c $(GEN_DIR)
	cd $(GEN_DIR) && $(MAKE)

FORCE:

clean: 
	rm -rf $(GEN_DIR) *.o driver

.PHONY: default clean

adrfunc.cpp:$(GEN_DIR)/libRapsodia.a
driver.cpp:$(GEN_DIR)/libRapsodia.a
