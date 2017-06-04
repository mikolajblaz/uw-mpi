# Autor: Mikołaj Błaż
# Nr indeksu: 346862

CC          := mpicc
CFLAGS      := -O3 -c -Wall -std=c99
LFLAGS      := -O3 -std=c99
OBJS        := collisions-common.o
HEADERS     := collisions-common.h
EXEC        := collisions-1 collisions-2 collisions-3
ALL         := $(EXEC)
TESTDIR     := myoff2

all : $(ALL)

collisions-1: collisions-1.o $(OBJS)
	$(CC) $(LFLAGS) -o $@ $^ -lm

collisions-2: collisions-2.o $(OBJS)
	$(CC) $(LFLAGS) -o $@ $^ -lm

collisions-3: collisions-3.o $(OBJS)
	$(CC) $(LFLAGS) -o $@ $^ -lm

%.o : %.c $(HEADERS)
	$(CC) $(CFLAGS) $<

clean :
	rm -f *.o *.out *.err $(ALL)

test: all
	mkdir -p tmp
	cp tests/$(TESTDIR)/gal* tmp
	cd tmp && \
	time mpirun -np 4 ../collisions-3 -v --hor 2 --ver 2 --gal1 gal1.txt --gal2 gal2.txt --delta 0.1 --total 1.0
	diff -r tmp tests/$(TESTDIR)
	echo -e "\n    OK!\n"
