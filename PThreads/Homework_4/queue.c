#include <stdio.h>
#include <stdlib.h>

#include "queue.h"

void Queue_construct(struct Queue *queue) {
    queue->first = NULL;
    queue->last = NULL;
    queue->size = 0;
}

void Queue_destruct(struct Queue *queue) {
    if (queue->size > 0) {
        struct QueueItem *current = queue->first;
        struct QueueItem *next = current->next;
        while (queue->size > 0) {
            free(current);
            current = next;
            next = current->next;
            queue->size--;
        }
    }
}

void *Queue_first(const struct Queue *queue) {
    if (queue->size == 0) {
        return NULL;
    }
    return queue->first->data;
}

void *Queue_last(const struct Queue *queue) {
    if (queue->size == 0) {
        return NULL;
    }
    return queue->last->data;
}

int Queue_pushFront(struct Queue *queue, void *data) {
    struct QueueItem *new_item =
        (struct QueueItem*) malloc(sizeof(struct QueueItem));

    if (new_item == NULL) {
        return 1; // Bad allocalion
    }

    new_item->prev = NULL;
    new_item->next = queue->first;
    new_item->data = data;

    if (queue->size == 0) {
        queue->last = new_item;
    } else /* queue->size > 0 */ {
        queue->first->prev = new_item;
    }
    queue->first = new_item;

    queue->size++;
    return 0;
}

int Queue_pushBack(struct Queue *queue, void *data) {
    struct QueueItem *new_item =
        (struct QueueItem*) malloc(sizeof(struct QueueItem));

    if (new_item == NULL) {
        return 1; // Bad allocalion
    }

    new_item->next = NULL;
    new_item->prev = queue->last;
    new_item->data = data;

    if (queue->size == 0) {
        queue->first = new_item;
    } else /* queue->size > 0 */ {
        queue->last->next = new_item;
    }
    queue->last = new_item;

    queue->size++;
    return 0;
}

int Queue_popFront(struct Queue *queue) {
    if (queue->size == 0) {
        return 1; // Empty queue
    }

    if (queue->size == 1) {
        free(queue->first);
        queue->first = NULL;
        queue->last = NULL;
    } else /* size > 1 */ {
        fflush(stdout);
        struct QueueItem *new_first = queue->first->next;
        new_first->prev = NULL;

        free(queue->first);
        queue->first = new_first;
    }
    queue->size--;
    return 0;
}

int Queue_popBack(struct Queue *queue) {
    if (queue->size == 0) {
        return 1; // Empty queue
    }

    if (queue->size == 1) {
        free(queue->last);
        queue->first = NULL;
        queue->last = NULL;
    } else /* size > 1 */ {
        struct QueueItem *new_last = queue->last->prev;
        new_last->next = NULL;

        free(queue->last);
        queue->last = new_last;
    }
    queue->size--;
    return 0;
}