CC=gcc
CFLAGS=-ofast-g-Wall
EXEC=projet
all:	$(EXEC)

$(EXEC):main.o simulation.o velocite.o	
	$(CC)	$^ -o $@ -lm

%.o:%.c
	$(CC) -c $(CFLAGS) $< -o $@

main.o:main.c simulation.h velocite.h	
lennardjones.o: simulation.c simulation.h
velocity.o: velocite.c	velocite.h
clean:
	rm -rf *.o $(EXEC)
