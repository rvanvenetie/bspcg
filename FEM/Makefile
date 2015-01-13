UNAME := $(shell whoami)
WARNING_FLAGS=-Wall -Wextra -Werror-implicit-function-declaration -Wshadow -Wstrict-prototypes -pedantic
CFLAGS= -g   -std=c99 -DMCBSP_COMPATIBILITY_MODE
IDIR = ../libs
LFLAGS=  -lpthread -lm  #-lrt

ifneq ($(UNAME),janner)
	CFLAGS += -pthread
endif

ifeq ($(UNAME),bissstud)
	CC = mpiicc
	LFLAGS += -lbsponmpi
else
	CC = gcc
	LFLAGS += $(IDIR)/compat-libmcbsp1.2.0.a 
endif


OBJ= fem.o $(IDIR)/mesh.o $(IDIR)/bspedupack.o $(IDIR)/io.o $(IDIR)/bspmv.o $(IDIR)/bspip.o $(IDIR)/vec.o
OBJ_BSP_FEM= bspfem.o  $(OBJ)

%.o: %.c
	$(CC) -I $(IDIR) -c $(CFLAGS) $<

all: bsp fem

fem: $(OBJ_BSP_FEM)
	$(CC) $(CFLAGS) -o bspfem $(OBJ_BSP_FEM) $(LFLAGS)

icrstest: $(OBJ) icrstest.o
	$(CC) $(CFLAGS) -o icrstest $(OBJ) icrstest.o $(LFLAGS)

bsp:
	cd $(IDIR); make all;

clean:
	rm -f *.o
	rm -f main
	rm -f m2mtx
	rm -f mtx2dm
	cd $(IDIR); make clean;
