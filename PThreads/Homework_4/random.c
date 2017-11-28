#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>

#include "random.h"
#include "log.h"

const size_t CHUNKS_COUNT = 5000;
const size_t CHUNK_SIZE = 5000;

void Generator_construct(struct Generator* generator) {
    generator->chunks = (int**) calloc(CHUNKS_COUNT, sizeof(int*));

    pthread_mutex_init(&generator->mutex, NULL);

    generator->state = RUN;
}

void Generator_destruct(struct Generator* generator) {
    for (int i = 0; i < CHUNKS_COUNT; ++i) {
        if (generator->chunks[i] != NULL) {
            free(generator->chunks[i]);
        }
    }
    free(generator->chunks);

    pthread_mutex_destroy(&generator->mutex);
}

void Chunk_construct(struct Chunk* chunk) {
    chunk->chunk = NULL;
    chunk->pos = 0;
}

void Chunk_destruct(struct Chunk* chunk) {
    if (chunk->chunk != NULL) {
       free(chunk->chunk);
    }
}

int* Generator_generateNewChunk() {
    int* chunk = (int*) calloc(CHUNK_SIZE, sizeof(int));
    for (int i = 0; i < CHUNK_SIZE; ++i) {
        chunk[i] = rand();
    }
    return chunk;
}

void Generator_run(struct Generator* generator) {
    int i = 0;
    while (generator->state == RUN) {
        while (generator->chunks[i] != NULL && generator->state == RUN) {
            sleep(1);
        }

        if (generator->state == STOP) {
            break;
        }
        generator->chunks[i] = Generator_generateNewChunk();

        i++;
        if (i >= CHUNKS_COUNT) {
            i = 0;
        }
    }
}

void Generator_stop(struct Generator* generator) {
    generator->state = STOP;
}

int* Generator_get_chunk(struct Generator* generator) {
    pthread_mutex_lock(&generator->mutex);

    while (generator->chunks[generator->pos] == NULL) {
    }

    int* result = generator->chunks[generator->pos];
    generator->chunks[generator->pos] = NULL;

    generator->pos++;
    if (generator->pos >= CHUNKS_COUNT) {
        generator->pos = 0;
    }

    pthread_mutex_unlock(&generator->mutex);

    return result;
}

int Generator_rand(struct Generator* generator, struct Chunk* chunk) {
    if (chunk->chunk == NULL) {
        chunk->chunk = Generator_get_chunk(generator);
        chunk->pos = 0;
    } else if (chunk->pos >= CHUNK_SIZE) {
        free(chunk->chunk);
        chunk->chunk = Generator_get_chunk(generator);
        chunk->pos = 0;
    }
    
    return chunk->chunk[chunk->pos++];
}