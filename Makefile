CC = mpicc
CFLAGS = -Wall -Wextra -O2 -fopenmp
LDFLAGS = -fopenmp

SOURCES = $(wildcard *.c)
HEADERS = $(wildcard *.h)
OBJ = $(SOURCES:.c=.o)

RUNFLAG = --hostfile $(OAR_NODEFILE)
C0 = --C0 54956c322f553ee3
C1 = --C1 ccb92d1a6ce6daec
N  = --n 15

all: main

run: main
	mpirun $(RUNFLAG) ./main $(N) $(C0) $(C1)

main: $(OBJ)
	$(CC) $(CFLAGS) $^ -o $@

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) main

.PHONY: all clean run
