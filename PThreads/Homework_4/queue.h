#ifndef QUEUE_H
#define QUEUE_H

#include <stdatomic.h>

struct QueueItem {
    void *next, *prev;
    void *data;
};

struct Queue {
    struct QueueItem *first, *last;
    atomic_size_t size;
};

void Queue_construct(struct Queue *queue);
void Queue_destruct(struct Queue *queue);

void *Queue_first(const struct Queue *queue);
void *Queue_last(const struct Queue *queue);

int Queue_pushFront(struct Queue *queue, void *data);
int Queue_pushBack(struct Queue *queue, void *data);

int Queue_popFront(struct Queue *queue);
int Queue_popBack(struct Queue *queue);

#endif // QUEUE_H