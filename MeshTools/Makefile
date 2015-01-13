UNAME := $(shell whoami)
WARNING_FLAGS=-Wall -Wextra -Werror-implicit-function-declaration -Wshadow -Wstrict-prototypes -pedantic
CFLAGS= -g   -std=c99 -DMCBSP_COMPATIBILITY_MODE -pthread
IDIR = ../libs
FDIR = ../FEM
LFLAGS=  -lpthread -lm  #-lrt

ifeq ($(UNAME),bissstud)
	CC = mpiicc
	LFLAGS += -lbsponmpi
else
	CC = gcc
	LFLAGS += $(IDIR)/compat-libmcbsp1.2.0.a 
endif

OBJ_M2FEM = m2fem.o $(IDIR)/mesh.o $(FDIR)/fem.o $(IDIR)/bspedupack.o $(IDIR)/io.o $(IDIR)/bspmv.o $(IDIR)/bspip.o $(IDIR)/vec.o
OBJ_M2MTX = m2mtx.o $(IDIR)/mesh.o
OBJ_MTX2DM =mtx2dm.o $(IDIR)/mesh.o

%.o: %.c
	$(CC) -I $(IDIR) -I $(FDIR) -c $(CFLAGS) $<

all:bsp fem m2mtx mtx2dm m2fem


m2fem: $(OBJ_M2FEM)
	$(CC) $(CFLAGS) -o m2fem $(OBJ_M2FEM) $(LFLAGS)

m2mtx: $(OBJ_M2MTX)
	$(CC) $(CFLAGS) -o m2mtx $(OBJ_M2MTX) $(LFLAGS)

mtx2dm: $(OBJ_MTX2DM)
	$(CC) $(CLFAGS) -o mtx2dm $(OBJ_MTX2DM) $(LFLAGS)

bsp:
	cd $(IDIR); make all;
fem:
	cd $(FDIR); make all;

clean:
	rm -f *.o
	rm -f m2fem
	rm -f m2mtx
	rm -f mtx2dm
	cd $(IDIR); make clean;
