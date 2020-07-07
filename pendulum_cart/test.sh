#!/bin/bash

# Compile & link
gcc -o run cartulum.c -O3 -I .. ../house.a \
	-lblas -llapacke -lsundials_cvode -lsundials_nvecserial -lm

# Run
./run

# Plot
gnuplot plot_cartulum.plt
