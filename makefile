CC=gcc
CFLAGS=-Wall -Wextra 
OFLAGS=-O2 -march=native -mtune=native
LFLAGS=-lm 

TARGET=main.exe

.PHONY: all clean

all: $(TARGET)

$(TARGET): main.o simulation.o
	$(CC) $(LFLAGS) $^ -o $@ -lm

%.o: %.c
	$(CC) -c $(CFLAGS) $(OFLAGS) $< -o $@

main.c: simulation.c simulation.h

matrix.c: simulation.h

clean:
	rm -Rf *~ *.o $(TARGET)
