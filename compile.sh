#!/bin/bash

gcc -c -O3 -Wall -g sigma.c

gcc -c -O3 -Wall -g statmat.c

gcc -c -O3 -Wall -g propagator.c

gcc -c -O3 -Wall -g pred_upd.c

gcc -c -O3 -Wall -g house_test.c

gcc -c -O3 -Wall -g noisemaker.c

ar rcs house.a *.o

#gcc -o test *.o -lblas -llapacke -lsundials_cvode -lsundials_nvecserial -lm

rm *.o
