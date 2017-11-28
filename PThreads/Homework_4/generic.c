#include <math.h>
#include <time.h>
#include <stdio.h>
#include <getopt.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <stdatomic.h>

#include "log.h"
#include "graph.h"
#include "random.h"
#include "threadpool.h"

#define SWAP(type, i, j) {type t = i; i = j; j = t;}
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

void get_time(struct timespec *to_store) {
    clock_gettime(CLOCK_MONOTONIC, to_store);
}

double elapsed_time(struct timespec from) {
    struct timespec now;
    get_time(&now);

    double elapsed = (now.tv_sec - from.tv_sec);
    elapsed += (now.tv_nsec - from.tv_nsec) / 1000000000.0;
    return elapsed;
}

enum Mode { NONE, GENERATE, FROM_FILE } mode = NONE;

struct ThreadPool thread_pool;
struct Generator generator;
struct Chunk** chunks;

struct graph_t* graph;

pthread_cond_t initialised_cond;
pthread_mutex_t initialised_mutex;

struct FitnessData {
    int (*fitness_func)(const int*, int, void*);
    void* data;
};

struct Genom {
    int* genom;
    int size;

    int fitness;
};

struct Population {
    struct Genom** genoms;
    int size;
    int genom_size;
};

void print_population(const struct Population* population) {
    printf("\n");
    for (int i = 0; i < population->size; ++i) {
        printf("%d\t:", i);
        for (int j = 0; j < population->genoms[i]->size; ++j) {
            printf(" %d", population->genoms[i]->genom[j]);
        }
        printf("\n");
    }
} 

void Genom_construct(struct Genom* genom, int size) {
    genom->size = size;
    genom->genom = (int*) calloc(genom->size, sizeof(int));
    genom->fitness = 0;
}

void Genom_destruct(struct Genom* genom) {
    free(genom->genom);
}

void Population_construct(struct Population* population, int population_size, int genom_size) {
    population->size = population_size;
    population->genom_size = genom_size;

    population->genoms = (struct Genom**) calloc(population_size, sizeof(struct Genom*));
    for (int i = 0; i < population_size; ++i) {
        population->genoms[i] = (struct Genom*) malloc(sizeof(struct Genom));
        Genom_construct(population->genoms[i], genom_size);
    }
}

void Population_destruct(struct Population* population) {
    for (int i = 0; i < population->size; i++) {
        Genom_destruct(population->genoms[i]);
        free(population->genoms[i]);
    }
    free(population->genoms);
}

struct ThreadRandom {
    struct Population* population;
    int population_size;

    const char* input_file;
    int generate_size;
};

int fitness_function(const int* genom, int size, void* data) {
    struct graph_t* graph = (struct graph_t*) data;
    int result = graph_weight(graph, genom[size - 1], genom[0]);

    for (int i = 1; i < size; ++i) {
        result += graph_weight(graph, genom[i - 1], genom[i]);
    }
    return result;
}

void* thread_random(void* arg) {
    struct ThreadRandom* params = (struct ThreadRandom*) arg;

    srand(time(NULL));

    if (mode == GENERATE) {
        graph = graph_generate(params->generate_size, 1000);
    } else {
        FILE* input_file = fopen(params->input_file, "r");
        graph = graph_read(input_file);
        fclose(input_file);
    }

    struct Population* population = params->population;
    Population_construct(params->population, params->population_size, graph->n);

    int* nums = (int*) calloc(graph->n, sizeof(int));
    for (int i = 0; i < params->population->size; ++i) {
        for (int j = 0; j < graph->n; ++j) {
            nums[j] = j;
        }
        for (int j = graph->n - 1; j >= 0; --j) {
            int r = rand() % (j + 1);
            population->genoms[i]->genom[j] = nums[r];

            SWAP(int, nums[r], nums[j]);
        }
        population->genoms[i]->fitness = fitness_function(population->genoms[i]->genom, 
                                                          population->genom_size, 
                                                          (void*) graph);
    }
    free(nums);

    Generator_construct(&generator);

    pthread_mutex_lock(&initialised_mutex);
    pthread_cond_signal(&initialised_cond);
    pthread_mutex_unlock(&initialised_mutex);

    Generator_run(&generator);
    Generator_destruct(&generator);

    free(params);
    return NULL;
}

struct ThreadSelect {
    const struct Population* population;
    struct Genom* new_genom;
};

void thread_select(void* arg, int thread_num) {
    struct ThreadSelect* params = (struct ThreadSelect*) arg;

    const struct Population* population = params->population;
    struct Genom* new_genom = params->new_genom;
    struct Chunk* chunk = chunks[thread_num];

    int pos_1 = Generator_rand(&generator, chunk) % population->size;
    int pos_2 = Generator_rand(&generator, chunk) % (population->size - 1);
    pos_2 += pos_1 + 1;
    pos_2 %= population->size;

    if (population->genoms[pos_1]->fitness < population->genoms[pos_2]->fitness) {
        memcpy(new_genom->genom, population->genoms[pos_1]->genom, population->genom_size * sizeof(int));
        new_genom->fitness = population->genoms[pos_1]->fitness; 
    } else if (population->genoms[pos_1]->fitness > population->genoms[pos_2]->fitness) {
        memcpy(new_genom->genom, population->genoms[pos_2]->genom, population->genom_size * sizeof(int));
        new_genom->fitness = population->genoms[pos_2]->fitness; 
    } else {
        int coin = Generator_rand(&generator, chunk) % population->size;
        if (coin == 0) {
            memcpy(new_genom->genom, population->genoms[pos_1]->genom, population->genom_size * sizeof(int));
            new_genom->fitness = population->genoms[pos_1]->fitness; 
        } else {
            memcpy(new_genom->genom, population->genoms[pos_2]->genom, population->genom_size * sizeof(int));
            new_genom->fitness = population->genoms[pos_2]->fitness; 
        }
    }

    free(params);
}

struct Population* selection(struct Population* population) {
    struct Population* result = (struct Population*) malloc(sizeof(struct Population));

    Population_construct(result, population->size, population->genom_size);

    for (int i = 0; i < population->size; ++i) {
        struct ThreadSelect* params = (struct ThreadSelect*) malloc(sizeof(struct ThreadSelect));
        params->population = population;
        params->new_genom = result->genoms[i];

        ThreadPool_add(&thread_pool, thread_select, (void*) params);
    }
    ThreadPool_wait(&thread_pool);

    Population_destruct(population);
    free(population);
    return result;
}

struct ThreadCross {
    struct Genom* parent_1;
    struct Genom* parent_2;
    struct Genom* child_1;
    struct Genom* child_2;
    struct FitnessData* data;
};

void thread_cross(void* arg, int thread_num) {
    struct ThreadCross* params = (struct ThreadCross*) arg;

    struct Genom* parent_1 = params->parent_1;
    struct Genom* parent_2 = params->parent_2;
    struct Genom* child_1 = params->child_1;
    struct Genom* child_2 = params->child_2;
    struct Chunk* chunk = chunks[thread_num];

    int genom_size = parent_1->size;

    int pos_1 = Generator_rand(&generator, chunk) % genom_size;
    int pos_2 = Generator_rand(&generator, chunk) % genom_size;

    int const_intervals[2][2], const_size;
    int move_intervals[2][2], move_size;

    if (pos_1 < pos_2) {
        const_intervals[0][0] = pos_1;
        const_intervals[0][1] = pos_2;
        const_size = 1;

        move_intervals[0][0] = 0;
        move_intervals[0][1] = pos_1;
        move_intervals[1][0] = pos_2;
        move_intervals[1][1] = genom_size;
        move_size = 2;
    } else if (pos_1 > pos_2) {
        move_intervals[0][0] = pos_2;
        move_intervals[0][1] = pos_1;
        move_size = 1;

        const_intervals[0][0] = 0;
        const_intervals[0][1] = pos_2;
        const_intervals[1][0] = pos_1;
        const_intervals[1][1] = genom_size;
        const_size = 2;
    } else /* pos_1 == pos_2 */ {
        move_size = 0;

        const_intervals[0][0] = 0;
        const_intervals[0][1] = genom_size;
        const_size = 1;
    }

    char* used = (char*) calloc(genom_size, sizeof(char));
    for (int num = 0; num < 2; ++num) {
        for (int i = 0; i < const_size; ++i) {
            for (int j = const_intervals[i][0]; j < const_intervals[i][1]; ++j) {
                child_1->genom[j] = parent_1->genom[j];
                used[child_1->genom[j]] = num + 1;
            }
        }

        int ind = 0;
        for (int i = 0; i < move_size; ++i) {
            for (int j = move_intervals[i][0]; j < move_intervals[i][1]; ++j) {
                while (used[parent_2->genom[ind]] == num + 1) {
                    ind++;
                }
                child_1->genom[j] = parent_2->genom[ind++];
                used[child_1->genom[j]] = num + 1;
            }
        }

        child_1->fitness = (*params->data->fitness_func)(child_1->genom, genom_size, params->data->data);

        SWAP(struct Genom*, parent_1, parent_2);
        SWAP(struct Genom*, child_1,  child_2);
    }
    free(used);

    free(params);
}

struct Population* crossover(struct Population* population, struct FitnessData* data) {
    struct Population* result = (struct Population*) malloc(sizeof(struct Population));

    Population_construct(result, population->size, population->genom_size);

    for (int i = 0; i < population->size - 1; i += 2) {

        struct ThreadCross* params = (struct ThreadCross*) malloc(sizeof(struct ThreadCross));

        params->parent_1 = population->genoms[i];
        params->parent_2 = population->genoms[i + 1];
        params->child_1 = result->genoms[i];
        params->child_2 = result->genoms[i + 1];
        params->data = data;

        ThreadPool_add(&thread_pool, thread_cross, (void*) params);
    }
    if (population->size & 1) {
        memcpy(result->genoms[result->size - 1]->genom, population->genoms[population->size - 1]->genom, population->genom_size * sizeof(int));
        result->genoms[result->size - 1]->fitness = population->genoms[population->size - 1]->fitness; 
    }
    ThreadPool_wait(&thread_pool);

    Population_destruct(population);
    free(population);
    return result;
}

struct ThreadMutate {
    struct Genom* genom;
};

void thread_mutate(void* arg, int thread_num) {
    struct ThreadMutate* params = (struct ThreadMutate*) arg;

    struct Genom* genom = params->genom;
    struct Chunk* chunk = chunks[thread_num];

    if (Generator_rand(&generator, chunk) % 100 == 0) {
        int pos_1 = Generator_rand(&generator, chunk) % genom->size;
        int pos_2 = Generator_rand(&generator, chunk) % (genom->size - 1);
        pos_2 += pos_1 + 1;
        pos_2 %= genom->size;

        SWAP(int, genom->genom[pos_1], genom->genom[pos_2]);
    }

    free(params);
}

void mutation(struct Population* population) {
    for (int i = 0; i < population->size; ++i) {
        struct ThreadMutate* params = (struct ThreadMutate*) malloc(sizeof(struct ThreadMutate));

        params->genom = population->genoms[i];

        ThreadPool_add(&thread_pool, thread_mutate, (void*) params);

    }
    ThreadPool_wait(&thread_pool);
}

struct Genom* find_best(struct Population* population) {
    int best_fitness = 0;
    for (int i = 1; i < population->size; ++i) {
        if (population->genoms[i]->fitness < population->genoms[best_fitness]->fitness) {
            best_fitness = i;
        }
    }
    return population->genoms[best_fitness];
}

struct Population* process(struct Population* population, int iteration_count, struct FitnessData* data) {
    for (int i = 0; i < iteration_count; ++i) {
        population = selection(population);
        population = crossover(population, data);
        mutation(population);
        // print_population(population);
    }
    return population;
}

void print_usage() {
    printf("usage: generic t n s [--generate NUM] [--file FILE]\n\nЗадача коммивояжера методом генетического алгоритма\n\noptional arguments:\n    -h              показать вспомогательное сообщение\n    --generate NUM  сгенерировать исходный граф автоматически\n    --file FILE     считать граф из заданого файла\n");
}

int main(int argc, char *argv[]){
    char* input_file = NULL;
    int generate_size = 0;

    int t, n, s;

    while (1) {
        int option_index = 0;
        static struct option long_options[] = {
            {"file",     required_argument, 0,  0},
            {"generate", required_argument, 0,  0},
            {0,          0,                 0,  0}
        };

        char c = getopt_long(argc, argv, "h", long_options, &option_index);
        if (c == -1) {
            break;
        }

        switch (c) {
            case 0:
                if (option_index == 0) {
                    input_file = optarg;
                    mode = FROM_FILE;
                } else if (option_index == 1) {
                    generate_size = atoi(optarg);
                    mode = GENERATE;
                }
                break;

            case 'h':
                print_usage();
                exit(EXIT_SUCCESS);
                break;

            default:
                print_usage();
                exit(EXIT_FAILURE);
        }
    }

    if (optind + 3 != argc) {
        printf("Wrong arguments\n");
        print_usage();
        exit(EXIT_FAILURE);
    }

    t = atoi(argv[optind++]);
    n = atoi(argv[optind++]);
    s = atoi(argv[optind++]);

    if (mode == NONE) {
        printf("No input graph\n");
        print_usage();
        exit(EXIT_FAILURE);
    } else if (mode == FROM_FILE && input_file == NULL) {
        exit(EXIT_FAILURE);
    } else if (mode == GENERATE && generate_size == 0) {
        exit(EXIT_FAILURE);
    }

    struct timespec start_time;
    get_time(&start_time);

    chunks = (struct Chunk**) calloc(t, sizeof(struct Chunk*));
    for (int i = 0; i < t; ++i) {
        chunks[i] = (struct Chunk*) malloc(sizeof(struct Chunk));
        Chunk_construct(chunks[i]);
    }

    struct Population* population = (struct Population*) malloc(sizeof(struct Population));
    struct ThreadRandom* params = (struct ThreadRandom*) malloc(sizeof(struct ThreadRandom));

    params->population = population;
    params->population_size = n;
    params->input_file = input_file;
    params->generate_size = generate_size;

    pthread_t random_thread;
    pthread_cond_init(&initialised_cond, NULL);
    pthread_mutex_init(&initialised_mutex, NULL);

    pthread_mutex_lock(&initialised_mutex);

    pthread_create(&random_thread, NULL, thread_random, (void*)params);

    pthread_cond_wait(&initialised_cond, &initialised_mutex);
    pthread_mutex_unlock(&initialised_mutex);

    ThreadPool_construct(&thread_pool, t);

    struct FitnessData* data = (struct FitnessData*) malloc(sizeof(struct FitnessData));
    data->fitness_func = fitness_function;
    data->data = (void*) graph;

    population = process(population, s, data);
    struct Genom* best_genom = find_best(population);
    int* way = (int*) calloc(best_genom->size, sizeof(int));
    memcpy(way, best_genom->genom, best_genom->size * sizeof(int));
    int cost = best_genom->fitness;

    free(data);

    ThreadPool_join(&thread_pool);
    ThreadPool_destruct(&thread_pool);

    Generator_stop(&generator);

    pthread_join(random_thread, NULL);
    
    pthread_mutex_destroy(&initialised_mutex);
    pthread_cond_destroy(&initialised_cond);
    
    graph_destroy(graph);
    Population_destruct(population);
    free(population);

    for (int i = 0; i < t; ++i) {
        Chunk_destruct(chunks[i]);
        free(chunks[i]);
    }
    free(chunks);

    double laps_time = elapsed_time(start_time);

    FILE* stats_file = fopen("stats.txt", "a");
    fprintf(stats_file, "%d %d %d %d %d %f %d\n", t, n, s, generate_size, s, laps_time, cost);
    for (int i = 0; i < generate_size; ++i) {
        fprintf(stats_file, "%d ", way[i]);
    }
    fprintf(stats_file, "\n");
    fclose(stats_file);

    free(way);

    exit(EXIT_SUCCESS);
    return 0;
}