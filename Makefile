
CC=g++ -std=c++11 -O2

all: exe_BN_driven exe_BNcc_driven exe_BN_STS

exe_BN_driven: BN_driven.cpp
	$(CC) -o exe_BN_driven BN_driven.cpp -lz

exe_BN_STS: BN_STS.cpp
	$(CC) -o exe_BN_STS BN_STS.cpp -lz

exe_BNcc_driven: BNcc_driven.cpp
	$(CC) -o exe_BNcc_driven BNcc_driven.cpp -lz

