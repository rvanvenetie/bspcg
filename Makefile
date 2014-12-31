CC= gcc
BSPDIR = ../../Linux-zoals-het-hoort/BSPedupack1.02/
WARNING_FLAGS=-Wall -Wextra -Werror-implicit-function-declaration -Wshadow -Wstrict-prototypes -pedantic
CFLAGS= -O3  -std=c99 -DMCBSP_COMPATIBILITY_MODE #-O3 -flibrary-level 2 -bspfifo 10000  -fcombine-puts -fcombine-puts-buffer 256K,128M,4K
LFLAGS= $(BSPDIR)/compat-libmcbsp1.2.0.a -pthread -lm -lrt
#Een of andere define hier voor jouw mac shit
IDIR = $(BSPDIR)

.c.o:
	$(CC) -I $(IDIR) -c $(CFLAGS) $<

all: main
main: main.o
	$(CC) $(CFLAGS) -o main main.o $(LFLAGS)


clean:
	rm -f *.o
	rm -f main

main.o     :   main.c
