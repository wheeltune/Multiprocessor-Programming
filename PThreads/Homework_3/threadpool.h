#ifndef THREADPOOL_H
#define THREADPOOL_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <string.h>
#include <stdatomic.h>

#include "queue.h"
#include "log.h"

struct ThreadPool {
    struct Queue tasks;

    pthread_mutex_t mutex;
    pthread_cond_t cond;

    pthread_t *threads;
    atomic_size_t num_threads;
    atomic_size_t num_waiting_threads;

    enum {
        WORKING, JOINING
    } mode;
};

struct ThreadPoolTask {
    void (*task)(void*);
    void *arg;
};

static void *ThreadPool_threadMain(void *arg);

void ThreadPool_construct(struct ThreadPool *thread_pool, size_t num_threads) {
    Queue_construct(&thread_pool->tasks);
    pthread_mutex_init(&thread_pool->mutex, NULL);
    pthread_cond_init(&thread_pool->cond, NULL);

    thread_pool->mode = WORKING;

    thread_pool->threads = (pthread_t*) calloc(sizeof(pthread_t), num_threads);
    thread_pool->num_threads = num_threads;
    thread_pool->num_waiting_threads = 0;
    for (int i = 0; i < num_threads; ++i) {
        pthread_create(&thread_pool->threads[i], NULL,
            ThreadPool_threadMain, (void*)thread_pool);
    }
}

void ThreadPool_destruct(struct ThreadPool *thread_pool) {
    free(thread_pool->threads);

    Queue_destruct(&thread_pool->tasks);
    pthread_cond_destroy(&thread_pool->cond);
    pthread_mutex_destroy(&thread_pool->mutex);
}

void ThreadPool_wait(struct ThreadPool *thread_pool) {
    while (thread_pool->num_waiting_threads < thread_pool->num_threads ||
            thread_pool->tasks.size > 0) {}
}

void ThreadPool_join(struct ThreadPool *thread_pool) {
    pthread_mutex_lock(&thread_pool->mutex);
    thread_pool->mode = JOINING;
    pthread_cond_broadcast(&thread_pool->cond);
    pthread_mutex_unlock(&thread_pool->mutex);

    for (int i = 0; i < thread_pool->num_threads; ++i) {
        pthread_join(thread_pool->threads[i], NULL);
    }
}

int ThreadPool_add(struct ThreadPool *thread_pool, void (*function)(void*), void *arg) {
    pthread_mutex_lock(&thread_pool->mutex);

    if (thread_pool->mode == JOINING) {
        return 1; // ThreadPool was shut downed
    }

    struct ThreadPoolTask *task =
        (struct ThreadPoolTask*) malloc(sizeof(struct ThreadPoolTask));
    task->task = function;
    task->arg = arg;

    Queue_pushBack(&thread_pool->tasks, task);
    pthread_cond_signal(&thread_pool->cond);

    pthread_mutex_unlock(&thread_pool->mutex);
    return 0;
}

static void *ThreadPool_threadMain(void *arg) {
    struct ThreadPool *thread_pool = (struct ThreadPool*) arg;

    for (;;) {
        pthread_mutex_lock(&thread_pool->mutex);

        while (thread_pool->mode == WORKING &&
               thread_pool->tasks.size == 0)
        {
            thread_pool->num_waiting_threads++;
            pthread_cond_wait(&thread_pool->cond, &thread_pool->mutex);
            thread_pool->num_waiting_threads--;
        }

        if (thread_pool->mode == JOINING && thread_pool->tasks.size == 0) {
            pthread_mutex_unlock(&thread_pool->mutex);
            break;
        }

        struct ThreadPoolTask *task = Queue_first(&thread_pool->tasks);
        Queue_popFront(&thread_pool->tasks);

        pthread_mutex_unlock(&thread_pool->mutex);

        ((*task->task))(task->arg);

        free(task);
    }

    return NULL;
}

#endif // THREADPOOL_H
