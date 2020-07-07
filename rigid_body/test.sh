#!/bin/bash

# Compile & link
gcc -o run rigid_body.c -O3 -I .. ../house.a \
        -lblas -llapacke -lsundials_cvode -lsundials_nvecserial -lm

# Run
./run

# Change directory
cd out

# Plot
gnuplot ../plot_rot.plt

# Generate PDF
pdflatex ../plots.tex

# Remove auxiliary files
rm *.aux
rm *.log

# Move PDF back to home directory
mv plots.pdf ../plots.pdf

# Back to home directory
cd ..
