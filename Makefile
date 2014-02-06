CC=g++-4.9
CFLAGS=--std=c++0x -O3 -fopenmp -DNDEBUG
INCLUDE=-I/usr/local/include
LINK=

all: gmres

gmres: gmres.cpp
	$(CC) $(CFLAGS) -o $@ $^ $(INCLUDE) $(LINK)

clean:
	$(RM) -f gmres_test *.o
