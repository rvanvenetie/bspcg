UNAME := $(shell whoami)
WARNING_FLAGS=-Wall -Wextra -Werror-implicit-function-declaration -Wshadow -Wstrict-prototypes -pedantic
CFLAGS= -g   -std=c99 -DMCBSP_COMPATIBILITY_MODE
IDIR = ../libs
LFLAGS=  -lpthread -lm  #-lrt

#ifeq ($(UNAME),raymond)
#	CC = mpicc
#	LFLAGS += -lbsponmpi
#else
	CC = gcc
	LFLAGS += $(IDIR)/compat-libmcbsp1.2.0.a 
#endif


OBJ= fem.o $(IDIR)/bspedupack.o $(IDIR)/io.o $(IDIR)/bspmv.o $(IDIR)/bspip.o $(IDIR)/vec.o

%.o: %.c
	$(CC) -I $(IDIR) -c $(CFLAGS) $<

all: bsp fem

fem: fem.o
	$(CC) $(CFLAGS) -o fem $(OBJ) $(LFLAGS)

bsp:
	cd $(IDIR); make all;

clean:
	rm -f *.o
	rm -f main
	cd $(IDIR); make clean;
