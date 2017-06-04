# Autor: Mikołaj Błaż
# Nr indeksu: 346862

CC          := mpicc
CFLAGS      := -O3 -c -Wall -std=c99
LFLAGS      := -O3 -std=c99
OBJS        := collisions.o collisions-helpers.o
HEADERS     := collisions-helpers.h collisions-optimizations.h
EXEC        := collisions-1 collisions-2 collisions-3
ALL         := $(EXEC)
TESTDIR     := myoff1

all : $(ALL)

collisions-1: collisions-optimizations1.o $(OBJS)
	$(CC) $(LFLAGS) -o $@ $^ -lm

collisions-2: collisions-optimizations2.o $(OBJS)
	$(CC) $(LFLAGS) -o $@ $^ -lm

collisions-3: collisions-optimizations3.o $(OBJS)
	$(CC) $(LFLAGS) -o $@ $^ -lm

collisions-optimizations1.o: collisions-optimizations.c $(HEADERS)
	$(CC) $(CFLAGS) -DOPTIMIZATION_1 -o collisions-optimizations1.o $<

collisions-optimizations2.o: collisions-optimizations.c $(HEADERS)
	$(CC) $(CFLAGS) -DOPTIMIZATION_2 -o collisions-optimizations2.o $<

collisions-optimizations3.o: collisions-optimizations.c $(HEADERS)
	$(CC) $(CFLAGS) -DOPTIMIZATION_3 -o collisions-optimizations3.o $<

%.o : %.c $(HEADERS)
	$(CC) $(CFLAGS) $<

clean :
	rm -f *.o *.out *.err $(ALL)

test: all
	mkdir -p tmp
	cp tests/$(TESTDIR)/gal* tmp
	cd tmp && \
	mpirun -np 10 ../collisions-1 -v --hor 5 --ver 2 --gal1 gal1.txt --gal2 gal2.txt --delta 0.1 --total 1.0
	diff -r tmp tests/$(TESTDIR)
	echo -e "\n    OK!\n"
	cd tmp && \
	mpirun -np 6 ../collisions-2 -v --hor 3 --ver 2 --gal1 gal1.txt --gal2 gal2.txt --delta 0.1 --total 1.0
	diff -r tmp tests/$(TESTDIR)
	echo -e "\n    OK!\n"
