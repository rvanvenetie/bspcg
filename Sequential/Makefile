CFLAGS = -Wall -g3 -I. -std=c99
LIBS = -framework Accelerate
OBJ = main.o
CC = gcc
default: $(OBJ)
	$(CC) -g -o main $^ $(LIBS) $(CFLAGS)

%.o: %.c $(DEPS)
	$(CC) -c -g -o $@ $< $(CFLAGS)

.PHONY: clean

clean:
	rm -f *.o main
	rm -rf *.dSYM
