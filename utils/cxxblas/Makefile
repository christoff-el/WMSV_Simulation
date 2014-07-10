OS = $(shell uname -s)

ifeq ($(OS),Darwin)
  #CC  := mpic++
  CC  := g++
endif

ifeq ($(OS),SunOS)
  CC  := g++
endif

ifeq ($(OS),Linux)
  CC  := g++
endif

CPPFLAGS = -framework vecLib -std=c++11 -Wall -O3 -I ~/proj/FLENS -c

.PHONY: all clean

all: ext_chol.o cxxblas_utils.a
ext_chol.o: ext_chol.cpp ext_chol.h
	$(CC) $(CPPFLAGS)  ext_chol.cpp
cxxblas_utils.a: ext_chol.o
	ar -rcu --target=mach-o-x86-64 libcxxblas_utils.a ext_chol.o 

clean:
	rm *.o libcxxblas_utils.a
