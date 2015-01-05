LFLAGS = compat-libmcbsp1.2.0.a -lm
CC = gcc
CFLAGS = -O3  -std=c99 -DMCBSP_COMPATIBILITY_MODE

io: io.o
	$(CC) $(CFLAGS) -o io io.o bspedupack.o $(LFLAGS)

io.o: bspedupack.o io.c io.h
	$(CC) $(CFLAGS) -c io.c

clean:
	rm -f *.o
	rm -f io

all: io
