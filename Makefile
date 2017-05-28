CC          := mpicc
CFLAGS      := -O3 -c -Wall
LFLAGS      := -O3
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
	mpirun -np 4 ./collisions -v --hor 2 --ver 2 --gal1 tests/in1gal1.txt --gal2 tests/in1gal2.txt --delta 0.1 --total 1.0
