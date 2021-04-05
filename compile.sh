g++ -c *.cpp --std=c++11 -Wall -pedantic -g -O3

rm -f house.a

ar rcs house.a *.o

rm -f *.o


