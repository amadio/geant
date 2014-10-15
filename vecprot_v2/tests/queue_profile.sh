#!/bin/bash
ARCH="gcc"`gcc -dumpversion`"("`gcc -dumpmachine`"_"`grep -c ^processor /proc/cpuinfo`"_core)"
echo "Compiling test_mpmc.cc using "$ARCH
g++ test_mpmc.cc -o test_mpmc -O3 -pthread -std=c++11 -I$HOME/boost/include -I`root-config --incdir` `root-config --libs`
rm -f q_*.txt
for NTHREADS in 1 2 3 4 5 6 7 8 12 16 24 32
do
  test_mpmc $NTHREADS $ARCH 5
done
root -b -q comparison.C+
