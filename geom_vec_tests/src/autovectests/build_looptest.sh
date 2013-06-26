g++ -ffast-math  -ftree-vectorize -ftree-vectorizer-verbose=6 -O2 -c t.cxx -o tgcc.o 2> logg++
icc -O2 -vec-report2 -c t.cxx -o ticc.o 2> logicc
