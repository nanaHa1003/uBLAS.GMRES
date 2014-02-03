CC=g++-4.9
CFLAGS=--std=c++0x -O3 -fopenmp
INCLUDE=-I/usr/local/include
LINK=-L/usr/local/lib

all: gmres_test

gmres_test: gmres.cpp
	$(CC) $(CFLAGS) -o $@ $^ $(INCLUDE) $(LINK)

clean:
	$(RM) -f gmres_test *.o
