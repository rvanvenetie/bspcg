UNAME := $(shell whoami)
WARNING_FLAGS=-Wall -Wextra -Werror-implicit-function-declaration -Wshadow -Wstrict-prototypes -pedantic
CFLAGS= -g   -std=c99 -DMCBSP_COMPATIBILITY_MODE
IDIR = ../libs
LFLAGS=  -pthread -lm  #-lrt

ifeq ($(UNAME),bissstud)
	CC = mpiicc
	LFLAGS += -lbsponmpi
else
	CC = gcc
	LFLAGS += $(IDIR)/compat-libmcbsp1.2.0.a 
endif


OBJ= mesh.o fem.o $(IDIR)/bspedupack.o $(IDIR)/io.o $(IDIR)/bspmv.o $(IDIR)/bspip.o $(IDIR)/vec.o

%.o: %.c
	$(CC) -I $(IDIR) -c $(CFLAGS) $<

all: bsp fem

fem: $(OBJ)
	$(CC) $(CFLAGS) -o fem $(OBJ) $(LFLAGS)

bsp:
	cd $(IDIR); make all;

mesh.o: mesh.c
	$(CC) -c $(CFLAGS) $<

m2mtx: mesh.o
	$(CC) $(CFLAGS) -o m2mtx mesh.o m2mtx.c

mtx2dm: mesh.o
	$(CC) $(CFLAGS) -o m2mtx mesh.o mtx2dm.c

clean:
	rm -f *.o
	rm -f main
	rm -f m2mtx
	rm -f mtx2dm
	cd $(IDIR); make clean;
