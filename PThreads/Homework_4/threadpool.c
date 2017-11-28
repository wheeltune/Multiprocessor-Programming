#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <string.h>
#include <unistd.h>
#include <stdatomic.h>

#include "queue.h"
#include "threadpool.h"

void ThreadPool_construct(struct ThreadPool *thread_pool, size_t num_threads) {
    Queue_construct(&thread_pool->tasks);

    pthread_mutex_init(&thread_pool->wait_mutex, NULL);
    pthread_cond_init(&thread_pool->wait_cond, NULL);

    pthread_mutex_init(&thread_pool->mutex, NULL);
    pthread_cond_init(&thread_pool->cond, NULL);

    thread_pool->mode = WORKING;

    thread_pool->threads = (pthread_t*) calloc(sizeof(pthread_t), num_threads);
    thread_pool->num_threads = num_threads;
    thread_pool->num_waiting_threads = 0;
    for (int i = 0; i < num_threads; ++i) {
        struct ThreadPool_Data* data = (struct ThreadPool_Data*) malloc(sizeof(struct ThreadPool_Data));
        data->thread_pool = thread_pool;
        data->thread_num = i;
        pthread_create(&thread_pool->threads[i], NULL,
            ThreadPool_thread_main, (void*)data);
    }
}

void ThreadPool_destruct(struct ThreadPool *thread_pool) {
    free(thread_pool->threads);

    Queue_destruct(&thread_pool->tasks);
    pthread_cond_destroy(&thread_pool->cond);
    pthread_mutex_destroy(&thread_pool->mutex);
}

void ThreadPool_wait(struct ThreadPool *thread_pool) {
    pthread_mutex_lock(&thread_pool->wait_mutex);
    while (thread_pool->num_waiting_threads < thread_pool->num_threads ||
            thread_pool->tasks.size > 0) {
        pthread_cond_wait(&thread_pool->wait_cond, &thread_pool->wait_mutex);
    }
    pthread_mutex_unlock(&thread_pool->wait_mutex);
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

int ThreadPool_add(struct ThreadPool *thread_pool, void (*function)(void*, int), void *arg) {
    pthread_mutex_lock(&thread_pool->mutex);

    if (thread_pool->mode == JOINING) {
        return 1; // thread_pool was shut downed
    }

    struct ThreadPool_Task *task =
        (struct ThreadPool_Task*) malloc(sizeof(struct ThreadPool_Task));
    task->task = function;
    task->params = arg;

    Queue_pushBack(&thread_pool->tasks, task);
    pthread_cond_signal(&thread_pool->cond);

    pthread_mutex_unlock(&thread_pool->mutex);
    return 0;
}

static void *ThreadPool_thread_main(void *arg) {
    struct ThreadPool_Data* thread_data = (struct ThreadPool_Data*) arg;
    struct ThreadPool* thread_pool = thread_data->thread_pool;

    for (;;) {
        pthread_mutex_lock(&thread_pool->mutex);

        while (thread_pool->mode == WORKING &&
               thread_pool->tasks.size == 0)
        {
            pthread_mutex_lock(&thread_pool->wait_mutex);
            thread_pool->num_waiting_threads++;
            pthread_cond_signal(&thread_pool->wait_cond);
            pthread_mutex_unlock(&thread_pool->wait_mutex);

            pthread_cond_wait(&thread_pool->cond, &thread_pool->mutex);
            
            pthread_mutex_lock(&thread_pool->wait_mutex);
            thread_pool->num_waiting_threads--;
            pthread_mutex_unlock(&thread_pool->wait_mutex);
        }

        if (thread_pool->mode == JOINING && thread_pool->tasks.size == 0) {
            pthread_mutex_unlock(&thread_pool->mutex);
            break;
        }

        struct ThreadPool_Task *task = Queue_first(&thread_pool->tasks);
        Queue_popFront(&thread_pool->tasks);

        pthread_mutex_unlock(&thread_pool->mutex);

        ((*task->task))(task->params, thread_data->thread_num);

        free(task);
    }

    free(thread_data);
    return NULL;
}