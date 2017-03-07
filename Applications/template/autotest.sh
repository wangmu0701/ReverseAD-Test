#!/bin/sh

#INSTANCE_PATH
export INSTANCE_PATH=./gmm
INPUTS=./../gmm/gmm_instances/1k/gmm_d20_K10.txt

#DEVORDER
export DEVORDER=2

#ADOLCMETHOD : 1:Direct, 2:Indirect, 3:FullHess, 4:Hess-V
export ADOLCMETHOD=3 

#REVERSEADMETHOD : 0:specific, 1:generic, 2:flat-code
export REVERSEADMETHOD=0

#LIBTAYLORMETHOD : 0:dense, 1:sparse
export LIBTAYLORMETHOD=0

#Compile executables
echo "Compiling executables..."
make -s -f Makefile.template

N=10
# run ADOL-C 10 times
for ((i=0;i<N;i++))
do
  ./test_adolc $INPUTS
done
for ((i=0;i<N;i++))
do
  ./test_reversead $INPUTS
done
for ((i=0;i<5;i++))
do
  ./test_libtaylor $INPUTS
done

make -s -f Makefile.template clean

unset INSTANCE_PATH
unset DEVORDER
unset ADOLCMETHOD
unset REVERSEADMETHOD
unset LIBTAYLORMETHOD

