#!/usr/bin/env bash

# ./DartThrowing
#   dimension
#   num_classes (positive for optimal rmatrix computation, negative for uniform off-diagonal entries)
#   priority (either c integers with low values indicating higher priority, or c floating points indicating class selection probability)
#   r_values (c*(c+1)/2 numbers in row major order of the upper matrix, or only c diagonal entries)
#   k_number (positive integer for the usual k number, negative integer for target number of samples, [0 1] float for rho-number, or positive float for specifying both the k-number/patience-factor and the rho-number)
#   domain_size (dimension numbers)

./DartThrowing 2 2 1 0 1 1 10 100 100
./DartThrowing 2 4 3 2 1 0 1 1 1 1 10 100 100
./DartThrowing 2 6 5 4 3 2 1 0 1 1 1 1 1 1 10 100 100
./DartThrowing 2 2 1 0 1 1 20 100 100
./DartThrowing 2 2 1 0 1 1 30 100 100

