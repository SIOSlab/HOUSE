rm -f out/*
rm -f plots/*

bash compile.sh

./ct.exe

./st.exe

gnuplot plotct.plt

gnuplot ploterr.plt
