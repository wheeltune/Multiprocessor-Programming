#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include "queue.h"
#include <stdatomic.h>

struct ThreadPool {
    struct Queue tasks;

    pthread_cond_t wait_cond;
    pthread_mutex_t wait_mutex;

    pthread_mutex_t mutex;
    pthread_cond_t cond;

    pthread_t* threads;
    atomic_size_t num_threads;
    atomic_size_t num_waiting_threads;

    enum {
        WORKING, JOINING
    } mode;
};

struct ThreadPool_Task {
    void (*task)(void*, int);
    void* params;
};

struct ThreadPool_Data {
    struct ThreadPool* thread_pool;
    int thread_num;
};

void ThreadPool_construct(struct ThreadPool *thread_pool, size_t num_threads);
void ThreadPool_destruct(struct ThreadPool *thread_pool);

void ThreadPool_wait(struct ThreadPool *thread_pool);
void ThreadPool_join(struct ThreadPool *thread_pool);
int  ThreadPool_add(struct ThreadPool *thread_pool, void (*function)(void*, int), void *arg);

static void *ThreadPool_thread_main(void *arg);

#endif // THREAD_POOL_H
