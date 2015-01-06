LFLAGS = compat-libmcbsp1.2.0.a -lm -lpthread
CC = gcc
CFLAGS = -O3  -std=c99 -DMCBSP_COMPATIBILITY_MODE

%.o : %.c
	$(CC) -c $(CFLAGS) $<

all: bspedupack.o io.o bspmv.o bspip.o vec.o


clean:
	rm -f *.o
	rm -f io

