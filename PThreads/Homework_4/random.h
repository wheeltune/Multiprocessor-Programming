#ifndef RANDOM_H
#define RANDOM_H

struct Generator {
    enum State {RUN, STOP} state;

    int** chunks;
    size_t pos;

    pthread_mutex_t mutex;
};

struct Chunk {
    int* chunk;
    size_t pos;
};

void Chunk_construct(struct Chunk* chunk);
void Chunk_destruct(struct Chunk* chunk);

void Generator_construct(struct Generator* generator);
void Generator_destruct(struct Generator* generator);

int* Generator_generateNewChunk();

void Generator_run(struct Generator* generator);
void Generator_stop(struct Generator* generator);

int* Generator_getChunk(struct Generator* generator);
int Generator_rand(struct Generator* generator, struct Chunk* chunk);

#endif // RANDOM_H