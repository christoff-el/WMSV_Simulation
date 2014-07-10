OS = $(shell uname -s)

ifeq ($(OS),Darwin)
  CC  := g++
endif

ifeq ($(OS),SunOS)
  CC  := g++
endif

ifeq ($(OS),Linux)
  CC  := g++
endif

CPPFLAGS = -DNDEBUG -DUSE_CXXLAPACK -framework veclib -std=c++11 -Wall -O3 -I ~/proj/FLENS -c

.PHONY: all clean

all: eye.o matmult.o matsqrt_diag.o trace.o pos_part.o norm1.o mat_solver.o expm.o inv.o mhg.o norminv.o det.o tracemul.o lib_utils.a

eye.o: eye.cpp
	$(CC) $(CPPFLAGS)  eye.cpp
	
matmult.o: matmult.cpp
	$(CC) $(CPPFLAGS)  matmult.cpp
	
matsqrt_diag.o: matsqrt_diag.cpp
	$(CC) $(CPPFLAGS)  matsqrt_diag.cpp
	
trace.o: trace.cpp
	$(CC) $(CPPFLAGS)  trace.cpp
pos_part.o: pos_part.cpp
	$(CC) $(CPPFLAGS)  pos_part.cpp
norm1.o: norm1.cpp
	$(CC) $(CPPFLAGS)  norm1.cpp
mat_solver.o: mat_solver.cpp
	$(CC) $(CPPFLAGS)  mat_solver.cpp
expm.o: expm.cpp
	$(CC) $(CPPFLAGS)  expm.cpp
inv.o: inv.cpp
	$(CC) $(CPPFLAGS)  inv.cpp
mhg.o: mhg.cpp
	$(CC) $(CPPFLAGS)  mhg.cpp
norminv.o: norminv.cpp
	$(CC) $(CPPFLAGS)  norminv.cpp
det.o: det.cpp
	$(CC) $(CPPFLAGS)  det.cpp
	
tracemul.o: tracemul.cpp
	$(CC) $(CPPFLAGS)  tracemul.cpp

	
lib_utils.a: eye.o matmult.o matsqrt_diag.o trace.o pos_part.o norm1.o mat_solver.o expm.o inv.o mhg.o norminv.o det.o tracemul.o
	ar -rcu --target=mach-o-x86-64 lib_utils.a eye.o matmult.o matsqrt_diag.o trace.o pos_part.o norm1.o mat_solver.o expm.o inv.o mhg.o norminv.o det.o tracemul.o

clean:
	rm *.o lib_utils.a
