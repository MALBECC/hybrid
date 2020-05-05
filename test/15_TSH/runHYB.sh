#!/bin/bash


file=ind7C
touch I7C.*; rm I7C.*; touch INPUT*; rm INPUT*
cp I7C_i.XV I7C.XV

export LIOHOME=path_to_lio_program
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LIOHOME/g2g:$LIOHOME/lioamber
export GFORTRAN_UNBUFFERED_ALL=1

HYB=../../bin
$HYB/hybrid < $file.fdf > $file.out

