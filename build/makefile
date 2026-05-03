CC      = gcc
CFLAGS  = -Wall -Wextra -pedantic -std=c89 -O2 -I../include

SRCDIR  = ../src
INCDIR  = ../include

all: pq_test_accuracy 

pq.o: $(SRCDIR)/pq.c $(INCDIR)/pq.h
	$(CC) $(CFLAGS) -c $(SRCDIR)/pq.c -o pq.o

pq_test_accuracy: $(SRCDIR)/pq_test_accuracy.c pq.o $(INCDIR)/pq.h
	$(CC) $(CFLAGS) $(SRCDIR)/pq_test_accuracy.c pq.o -o pq_test_accuracy -lm

test: test_accuracy 
	./pq_test_accuracy

clean:
	rm -f *.o pq_test_accuracy 