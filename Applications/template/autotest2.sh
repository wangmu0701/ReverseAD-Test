#!/bin/sh

#INSTANCE_PATH
export INSTANCE_PATH=./gmm

#DEVORDER
export DEVORDER=2

#ADOLCMETHOD : 1:Direct, 2:Indirect, 3:FullHess, 4:Hess-V, 5:DenseHess
export ADOLCMETHOD=6

#REVERSEADMETHOD : 0:specific, 1:generic, 2:flat-code
export REVERSEADMETHOD=2

#LIBTAYLORMETHOD : 0:dense, 1:sparse
export LIBTAYLORMETHOD=0

#Compile executables
echo "Compiling executables..."
make -s -f Makefile.template

N=10
INPUTS=./../gmm/gmm_instances/1k/gmm_d10_K5.txt
echo $INPUTS
for ((i=0;i<N;i++))
do
  ./test_adolc $INPUTS
done
INPUTS=./../gmm/gmm_instances/1k/gmm_d10_K10.txt
echo $INPUTS
for ((i=0;i<N;i++))
do
  ./test_adolc $INPUTS
done
INPUTS=./../gmm/gmm_instances/1k/gmm_d20_K5.txt
echo $INPUTS
for ((i=0;i<N;i++))
do
  ./test_adolc $INPUTS
done
INPUTS=./../gmm/gmm_instances/1k/gmm_d10_K25.txt
echo $INPUTS
for ((i=0;i<N;i++))
do
  ./test_adolc $INPUTS
done
INPUTS=./../gmm/gmm_instances/1k/gmm_d20_K10.txt
echo $INPUTS
for ((i=0;i<N;i++))
do
  ./test_adolc $INPUTS
done

make -s -f Makefile.template clean

unset INSTANCE_PATH
unset DEVORDER
unset ADOLCMETHOD
unset REVERSEADMETHOD
unset LIBTAYLORMETHOD

