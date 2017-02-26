#!/bin/sh

#INSTANCE_PATH
export INSTANCE_PATH=./mesh
INPUTS=./../mesh/gear.mesh

#DEVORDER
export DEVORDER=2

#ADOLCMETHOD : 1:Direct, 2:Indirect, 3:FullHess, 4:Hess-V
export ADOLCMETHOD=1 

#REVERSEADMETHOD : 0:specific, 1:generic, 2:flat-code
export REVERSEADMETHOD=0

#LIBTAYLORMETHOD : 0:dense, 1:sparse
export LIBTAYLORMETHOD=1

#Compile executables
echo "Compiling executables..."
make -s -f Makefile.template

N=1
# run ADOL-C 10 times
for ((i=0;i<N;i++))
do
  ./test_adolc $INPUTS
done
for ((i=0;i<N;i++))
do
  ./test_reversead $INPUTS
done
for ((i=0;i<N;i++))
do
  ./test_libtaylor $INPUTS
done

make -s -f Makefile.template clean

unset INSTANCE_PATH
unset DEVORDER
unset ADOLCMETHOD
unset REVERSEADMETHOD
unset LIBTAYLORMETHOD

