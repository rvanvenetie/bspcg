CC = gcc
WARNING_FLAGS=-Wall -Wextra -Werror-implicit-function-declaration -Wshadow -Wstrict-prototypes -pedantic
CFLAGS= -O3  -std=c99 -DMCBSP_COMPATIBILITY_MODE
IDIR = libs
LFLAGS= $(IDIR)/compat-libmcbsp1.2.0.a -lpthread -lm #-lrt

OBJ= main.o $(IDIR)/bspedupack.o $(IDIR)/io.o $(IDIR)/bspmv.o

%.o: %.c
	$(CC) -I $(IDIR) -c $(CFLAGS) $<

all: main

main: main.o
	$(CC) $(CFLAGS) -o main $(OBJ) $(LFLAGS)

bsp:
	cd libs; make all;

clean:
	rm -f *.o
	rm -f main
