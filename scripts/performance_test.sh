#!/usr/bin/env bash

# NOTE THAT THIS SCRIPT IS NOT FOR DIRECT EXECUTION!!!

# ./DartThrowing
#   dimension
#   num_classes (positive for optimal rmatrix computation, negative for uniform off-diagonal entries)
#   priority (either c integers with low values indicating higher priority, or c floating points indicating class selection probability)
#   r_values (c*(c+1)/2 numbers in row major order of the upper matrix, or only c diagonal entries)
#   k_number (positive integer for the usual k number, negative integer for target number of samples, [0 1] float for rho-number, or positive float for specifying both the k-number/patience-factor and the rho-number)
#   domain_size (dimension numbers)

# Performance test
./DartThrowing 2 2 1 0 1 1 10 100 100
./DartThrowing 2 4 3 2 1 0 1 1 1 1 10 100 100
./DartThrowing 2 6 5 4 3 2 1 0 1 1 1 1 1 1 10 100 100
./DartThrowing 2 2 1 0 1 1 20 100 100
./DartThrowing 2 2 1 0 1 1 30 100 100

# Correctness test & visualization
./DartThrowing 2 2 1 0 0.01 0.01 10 1 1 > ../testdata/pointset1.txt
./Main

# Demo for original algorithm
./DartThrowing 2 2 1 0 0.01 0.01 10 1 1 > ../testdata/pointset1.txt
./SFT ../testdata/pointset1.txt ../testdata/spectrum.txt 2 -1 1 1 256 -1
./SFT ../testdata/pointset1.txt ../testdata/spectrum0.txt 2 0 1 1 256 -1
./SFT ../testdata/pointset1.txt ../testdata/spectrum1.txt 2 1 1 1 256 -1
./PFM2PPM ../testdata/spectrum.txt ../testdata/spectrum.txt 1 1 0
./PFM2PPM ../testdata/spectrum0.txt ../testdata/spectrum0.txt 1 1 0
./PFM2PPM ../testdata/spectrum1.txt ../testdata/spectrum1.txt 1 1 0
./Main

# Demo for our algorithm
./ExDartThrowing 2 2 1 0 0.01 0.01 10 1 1 > ../testdata/pointset1.txt
./SFT ../testdata/pointset1.txt ../testdata/spectrum.txt 2 -1 1 1 256 -1
./SFT ../testdata/pointset1.txt ../testdata/spectrum0.txt 2 0 1 1 256 -1
./SFT ../testdata/pointset1.txt ../testdata/spectrum1.txt 2 1 1 1 256 -1
./PFM2PPM ../testdata/spectrum.txt ../testdata/spectrum.txt 1 1 0
./PFM2PPM ../testdata/spectrum0.txt ../testdata/spectrum0.txt 1 1 0
./PFM2PPM ../testdata/spectrum1.txt ../testdata/spectrum1.txt 1 1 0
./Main
