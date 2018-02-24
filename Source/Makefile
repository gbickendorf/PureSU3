CC=g++ 
CXXFLAGS = -W -Wall -Wextra -Wsign-compare \
					 --std=c++11 \
					 -m64 \
					 -msse \
					 -O3 \
					 -march=native \
					 -funroll-loops \
					 -ftree-vectorize \
					 -fomit-frame-pointer\
					 -fopenmp
main :	TwoLinkOperator.o\
		Matrix.o\
		SUNLattice.o\
		Random.o\
		Vectors.o\
		LatticeSettings.o
clean:
	rm -f *.o
