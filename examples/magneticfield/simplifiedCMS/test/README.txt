Test Benchmarking.cpp.AvgOverRuns : Averages time over 1/16 portion of total runs. Reports mean time and standard deviation.

CheckVectorAgainstScalar.cpp : Scalar version is verified for correctness using MainMagFieldTest.cpp. This file checks the correctness of vector version using scalar values as reference. If assert fails anywhere after k>16, then vector version is right. Basically verified that indexing in vector version is right. The failure of assert is due to very small relative difference later. Reason for this : yet unknown.

FloatDoubleDiff.cpp : Uses double template MagField.h to find relative error when floats are used instead of doubles. Benchmarked against values obtained using double. Tells about loss of precision.

