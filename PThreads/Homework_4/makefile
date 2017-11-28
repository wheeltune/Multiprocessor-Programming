OBJS = generic.o queue.o threadpool.o log.o graph.o random.o
CC = gcc
DEBUG = -g -pthread
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

run : $(OBJS) generic.c
	$(CC) $(LFLAGS) $(OBJS) -o generic

generic.o : generic.c threadpool.h log.h random.h
	$(CC) $(CFLAGS) generic.c

treadpool.o : threadpool.c threadpool.h queue.h
	$(CC) $(CFLAGS) threadpool.c

random.o : random.c random.h log.h
	$(CC) $(CFLAGS) random.c

queue.o : queue.h queue.c
	$(CC) $(CFLAGS) queue.c

log.o : log.h log.c
	$(CC) $(CFLAGS) log.c

graph.o : graph.h graph.c
	$(CC) $(CFLAGS) graph.c

clean:
	rm -r *.o *.dSYM generic
