# /***************************** Makefile ******************************/

SOURCES = $(wildcard *.c)
HEADERS = $(wildcard *.h)
OBJ = $(SOURCES:.c=.o)

FLAGS = -Wall -fopenmp

all: main

main: $(OBJ)
	gcc $(FLAGS) $^ -o $@

%.o: %.c $(HEADERS)
	gcc $(FLAGS) -c $< -o $@

clean:
	rm -f *.o main

.PHONY: all clean
