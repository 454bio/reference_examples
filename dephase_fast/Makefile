CC=g++
CXXFLAGS = -std=c++11 -O3
CFLAGS=-I.
DEPS = Model4.h
OBJ = Model4.o Basecall.o 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

Basecall: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

