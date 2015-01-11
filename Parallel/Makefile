UNAME := $(shell whoami)
WARNING_FLAGS=-Wall -Wextra -Werror-implicit-function-declaration -Wshadow -Wstrict-prototypes -pedantic
CFLAGS= -g   -std=c99 -DMCBSP_COMPATIBILITY_MODE
IDIR = ../libs
LFLAGS=  -lpthread -lm  #-lrt

ifeq ($(UNAME),raymond)
	CC = mpicc
	LFLAGS += -lbsponmpi
else
	CC = gcc
	LFLAGS += $(IDIR)/compat-libmcbsp1.2.0.a 
endif


OBJ= main.o $(IDIR)/bspedupack.o $(IDIR)/io.o $(IDIR)/bspmv.o $(IDIR)/bspip.o $(IDIR)/vec.o

%.o: %.c
	$(CC) -I $(IDIR) -c $(CFLAGS) $<

all: bsp main

main: main.o
	$(CC) $(CFLAGS) -o main $(OBJ) $(LFLAGS)

gen: 
	g++ rand_sym_mat.cpp -o rand_sym_mat
bsp:
	cd $(IDIR); make all;

clean:
	rm -f *.o
	rm -f main
	cd $(IDIR); make clean;
