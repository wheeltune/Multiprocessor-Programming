#include <pthread.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include "threadpool.c"

void test(void *arg) {
    sleep(1);
    printf("%ld HELLO WORLD\n", (long) arg);
}

int main(int argc, char* argv[]) {
    struct ThreadPool thread_pool;
    ThreadPool_construct(&thread_pool, 4);

    for (long i = 0; i < 10; ++i) {
        ThreadPool_add(&thread_pool, test, (void*) i);
        if (i % 2 == 1) {   
            ThreadPool_wait(&thread_pool);
        }
    }

    ThreadPool_join(&thread_pool);
    ThreadPool_destruct(&thread_pool);

    return 0;
}
