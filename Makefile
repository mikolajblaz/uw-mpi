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
