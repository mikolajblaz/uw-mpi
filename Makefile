CC          := mpicc
CFLAGS      := -O3 -c -Wall -std=c99
LFLAGS      := -O3 -std=c99
OBJS        := collisions.o collisions-helpers.o
HEADERS     := collisions-helpers.h
EXEC        := collisions
ALL         := $(EXEC)

all : $(ALL)

$(EXEC) : $(OBJS)
	$(CC) $(LFLAGS) -o $@ $^ -lm

%.o : %.c $(HEADERS)
	$(CC) $(CFLAGS) $<

clean :
	rm -f *.o *.out *.err $(ALL)

test: all
	mpicc -O3 -c -Wall -DDEBUG collisions.c
	mpicc -O3 -o collisions collisions.o collisions-helpers.o -lm
	cp tests/smoke/gal* tmp
	cd tmp && \
	mpirun -np 4 ../collisions -v --hor 2 --ver 2 --gal1 gal1.txt --gal2 gal2.txt --delta 0.1 --total 1.0
	diff -r tmp tests/smoke
