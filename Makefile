CC = mpicc
CFLAGS = -Wall -fopenmp -O2
LDFLAGS = -fopenmp
SOURCES = $(wildcard *.c)
HEADERS = $(wildcard *.h)
OBJ = $(SOURCES:.c=.o)

all: main

main: $(OBJ)
	$(CC) $(LDFLAGS) $^ -o $@

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) main

.PHONY: all clean
