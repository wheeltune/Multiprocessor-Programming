#include <mpi.h>
#include <omp.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#define MASTER 0
#define CHUNK_LENGTH 500
#define CHUNKS_COUNT 500

enum side_e {
    LEFT = 0, RIGHT = 1, UP = 2, DOWN = 3, NONE = 4
};

struct point_t {
    int x, y;
    int r;
};

int point_t_cmp (const void* a, const void* b) {
    struct point_t* p1 = (struct point_t*) a;
    struct point_t* p2 = (struct point_t*) b;
    if (p1->y == p2->y) {
        if (p1->x == p2->x) {
            return p1->r - p2->r;
        }
        return p1->x - p2->x;
    }
    return p1->y - p2->y;
}

struct chunk_t {
    int* data;
    int pos;
};

void chunk_init(struct chunk_t* chunk) {
    chunk->data = malloc(CHUNK_LENGTH * sizeof(int));
    chunk->pos = CHUNK_LENGTH;
}

void chunk_destroy(struct chunk_t* chunk) {
    free(chunk->data);
}

int l_, a_, b_, N_;

double start_time;
void run_generator(int count_fields) {
    MPI_File data_file;
    MPI_File_open(MPI_COMM_WORLD, "data.bin", MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &data_file);
    MPI_File_close(&data_file);

    int** chunks_data = (int**) calloc(CHUNKS_COUNT, sizeof(int*));
    int count_finished = 0;
    #pragma omp parallel shared(count_finished, chunks_data) num_threads(2)
    {
        if (omp_get_thread_num() == 0) {
            srand(time(NULL));
            int* finished = (int*) calloc(count_fields, sizeof(int));
            int chunk_pos = 0;

            int interlocutor_rank = 1;
            while (count_finished < count_fields) {
                if (!finished[interlocutor_rank - 1]) {
                    int flag = 0;
                    MPI_Status status;
                    MPI_Iprobe(interlocutor_rank, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
                    if (flag) {
                        int is_finished = 0;
                        MPI_Recv(&is_finished, 1, MPI_INT, interlocutor_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        if (is_finished) {
                            #pragma omp atomic
                            count_finished += 1;

                            finished[interlocutor_rank - 1] = 1;
                            continue;
                        }

                        while (chunks_data[chunk_pos] == NULL) {}
                        MPI_Ssend(chunks_data[chunk_pos], CHUNK_LENGTH, MPI_INT, interlocutor_rank, 0, MPI_COMM_WORLD);
                        free(chunks_data[chunk_pos]);
                        chunks_data[chunk_pos] = NULL;

                        chunk_pos++;
                        if (chunk_pos == CHUNKS_COUNT) {
                            chunk_pos = 0;
                        }
                    }
                }
                interlocutor_rank++;
                if (interlocutor_rank > count_fields) {
                    interlocutor_rank = 1;
                }
            }

            free(finished);

         } else {
            int chunk_pos = 0;
            while (count_finished < count_fields) {
                int* new_chunk = (int*) malloc(CHUNK_LENGTH * sizeof(int));
                for (int i = 0; i < CHUNK_LENGTH; ++i) {
                    new_chunk[i] = rand();
                }

                while (chunks_data[chunk_pos] != NULL && count_finished < count_fields) {}
                chunks_data[chunk_pos] = new_chunk;

                chunk_pos++;
                if (chunk_pos >= CHUNKS_COUNT) {
                    chunk_pos = 0;
                }
            }
        }
    }

    for (int i = 0; i < CHUNKS_COUNT; ++i) {
        if (chunks_data[i] != NULL) {
            free(chunks_data[i]);
        }
    }
    free(chunks_data);
    double laps_time = MPI_Wtime() - start_time;

    FILE* stats_file = fopen("stats.txt", "w");
    fprintf(stats_file, "%d %d %d %d %lfs\n", l_, a_, b_, N_, laps_time);
    fclose(stats_file);
}

int get_random(struct chunk_t* chunk) {
    if (chunk->pos >= CHUNK_LENGTH) {
        int is_finished = 0;
        MPI_Ssend(&is_finished, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);

        MPI_Status status;
        MPI_Recv(chunk->data, CHUNK_LENGTH, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        chunk->pos = 0;
    }
    return chunk->data[chunk->pos++];
}

double get_random_double(struct chunk_t* chunk) {
    return (double) get_random(chunk) / (double) RAND_MAX;
}

MPI_Offset get_offset(struct point_t* p) {
    return (p->y * l_ * a_ * a_ * b_ + p->x * a_ * b_  + p->r) * sizeof(int);
}

void run_process(int rank, int count_fields) {
    MPI_File data_file;
    MPI_File_open(MPI_COMM_WORLD, "data.bin", MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &data_file);

    struct point_t min_point, max_point;
    int horizontal_index = (rank - 1) % a_;
    int vertical_index = (rank - 1) / a_;
    min_point.x = horizontal_index * l_;
    min_point.y = vertical_index * l_;
    min_point.r = 0;
    max_point.x = min_point.x + l_ - 1;
    max_point.y = max_point.y + l_ - 1;
    max_point.r = count_fields;

    struct chunk_t chunk;
    chunk_init(&chunk);
    struct point_t* points = (struct point_t*) malloc(N_ * sizeof(struct point_t));
    int count_points = N_;
    for (int i = 0; i < N_; ++i) {
        points[i].x = get_random(&chunk) % l_ + min_point.x;
        points[i].y = get_random(&chunk) % l_ + min_point.y;
        points[i].r = get_random(&chunk) % count_fields;
    }
    qsort(points, N_, sizeof(struct point_t), point_t_cmp);

    int i = 0;
    while (i < N_) {
        int j = i;
        while (j < N_ && points[j].r == points[i].r) {
            j++;
        }
        int count = j - i;
        MPI_File_write_at(data_file, get_offset(points + i), &count, 1, MPI_INT, MPI_STATUS_IGNORE);

        i = j;
    }

    MPI_File_close(&data_file);

    {
        int is_finished = 1;
        MPI_Ssend(&is_finished, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
    }
    free(points);
    chunk_destroy(&chunk);
}

int main(int argc, char **argv) {
    l_ = atoi(argv[1]);
    a_ = atoi(argv[2]), b_ = atoi(argv[3]);
    N_ = atoi(argv[4]);

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    start_time = MPI_Wtime();
    if (rank == MASTER) {
        run_generator(size - 1);
    } else {
        run_process(rank, size - 1);
    }

    MPI_Finalize();
    return 0;
}
