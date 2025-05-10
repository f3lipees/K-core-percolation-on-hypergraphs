/**************************************************************************************************
 * If you use this code, please cite:
 * G. Bianconi and S. N. Dorogovstev
 * "Nature of hypergraph k-core percolation problems"
 * Physical Review E, 109, 014307 (2024).
 *
 * Version: 2.1
 *
 **************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <errno.h>
#include <limits.h>

void init_random(unsigned long seed);
double random_double(void);

typedef struct {
    int num_nodes;              // Number of nodes
    int num_hyperedges;         // Number of hyperedges
    int k_core;                 // k of the core
    int s_core;                 // s of the core (called Ncore in original)
    float avg_membership;       // Average membership (mav)
    int num_simulations;        // Number of Monte Carlo simulations
    float error_tolerance;      // Error tolerance
    char output_filename[256];  // Output filename
    unsigned long random_seed;  // Random seed (0 means use time as seed)
    bool verbose;               // Enable verbose output
} config_t;

typedef struct {
    int num_nodes;
    int num_hyperedges;
    int total_elements;         

    int* node_degrees;          
    int* hyperedge_sizes;       
    int** adjacency_lists;      
    float** weights;           

    int* visited;
    int* damaged;
    int* temp_state;
    int* cluster_sizes;
} hypergraph_t;

bool parse_arguments(int argc, char** argv, config_t* config);
void print_help();
hypergraph_t* create_hypergraph(const config_t* config);
void free_hypergraph(hypergraph_t* hg);
bool generate_random_hypergraph(hypergraph_t* hg, const config_t* config);
int find_connected_component(hypergraph_t* hg, int start_node, int cluster_id, const config_t* config);
bool run_percolation_simulation(hypergraph_t* hg, const config_t* config, float* results);
bool write_results_to_file(const float* results, int num_points, const config_t* config);
void print_progress(int current, int total);

// Platform-independent random number utilities
// To replace drand48 and srand48 which are not available on Windows
static unsigned long next_random = 1;

void init_random(unsigned long seed) {
    next_random = seed ? seed : (unsigned long)time(NULL);
}

double random_double(void) {
    // Simple linear congruential generator
    next_random = next_random * 1103515245 + 12345;
    return ((double)((next_random / 65536) % 32768)) / 32767.0;
}

void set_default_config(config_t* config) {
    config->num_nodes = 3000;
    config->num_hyperedges = 6000;
    config->k_core = 2;
    config->s_core = 2;
    config->avg_membership = 2.5;
    config->num_simulations = 100;
    config->error_tolerance = 0.01;
    strncpy(config->output_filename, "hypergraph_kcore_results.txt", sizeof(config->output_filename) - 1);
    config->random_seed = 0;  
    config->verbose = false;
}

bool parse_arguments(int argc, char** argv, config_t* config) {
    set_default_config(config);

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_help();
            return false;
        } else if (strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "--nodes") == 0) {
            if (++i < argc) {
                int val = atoi(argv[i]);
                if (val > 0) {
                    config->num_nodes = val;
                } else {
                    fprintf(stderr, "Error: Number of nodes must be positive\n");
                    return false;
                }
            }
        } else if (strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "--hyperedges") == 0) {
            if (++i < argc) {
                int val = atoi(argv[i]);
                if (val > 0) {
                    config->num_hyperedges = val;
                } else {
                    fprintf(stderr, "Error: Number of hyperedges must be positive\n");
                    return false;
                }
            }
        } else if (strcmp(argv[i], "-k") == 0 || strcmp(argv[i], "--kcore") == 0) {
            if (++i < argc) {
                int val = atoi(argv[i]);
                if (val > 0) {
                    config->k_core = val;
                } else {
                    fprintf(stderr, "Error: k-core parameter must be positive\n");
                    return false;
                }
            }
        } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--score") == 0) {
            if (++i < argc) {
                int val = atoi(argv[i]);
                if (val > 0) {
                    config->s_core = val;
                } else {
                    fprintf(stderr, "Error: s-core parameter must be positive\n");
                    return false;
                }
            }
        } else if (strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--membership") == 0) {
            if (++i < argc) {
                float val = atof(argv[i]);
                if (val > 0) {
                    config->avg_membership = val;
                } else {
                    fprintf(stderr, "Error: Average membership must be positive\n");
                    return false;
                }
            }
        } else if (strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "--runs") == 0) {
            if (++i < argc) {
                int val = atoi(argv[i]);
                if (val > 0) {
                    config->num_simulations = val;
                } else {
                    fprintf(stderr, "Error: Number of simulations must be positive\n");
                    return false;
                }
            }
        } else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) {
            if (++i < argc) {
                strncpy(config->output_filename, argv[i], sizeof(config->output_filename) - 1);
            }
        } else if (strcmp(argv[i], "--seed") == 0) {
            if (++i < argc) {
                unsigned long val = strtoul(argv[i], NULL, 10);
                config->random_seed = val;
            }
        } else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            config->verbose = true;
        } else {
            fprintf(stderr, "Error: Unknown option: %s\n", argv[i]);
            print_help();
            return false;
        }
    }

    return true;
}

void print_help() {
    printf("Hypergraph K-Core Percolation Simulator\n\n");
    printf("Usage: hypergraph_kcore [options]\n\n");
    printf("Options:\n");
    printf("  -h, --help              Display this help message\n");
    printf("  -n, --nodes N           Set number of nodes (default: 3000)\n");
    printf("  -e, --hyperedges E      Set number of hyperedges (default: 6000)\n");
    printf("  -k, --kcore K           Set k-core parameter (default: 2)\n");
    printf("  -s, --score S           Set s-core parameter (default: 2)\n");
    printf("  -m, --membership M      Set average membership (default: 2.5)\n");
    printf("  -r, --runs R            Set number of Monte Carlo simulations (default: 100)\n");
    printf("  -o, --output FILE       Set output filename (default: hypergraph_kcore_results.txt)\n");
    printf("      --seed SEED         Set random seed (default: use current time)\n");
    printf("  -v, --verbose           Enable verbose output\n");
}

hypergraph_t* create_hypergraph(const config_t* config) {
    hypergraph_t* hg = NULL;

    // Allocate the main structure
    hg = (hypergraph_t*)calloc(1, sizeof(hypergraph_t));
    if (!hg) {
        fprintf(stderr, "Error: Memory allocation failed for hypergraph structure\n");
        return NULL;
    }

    hg->num_nodes = config->num_nodes;
    hg->num_hyperedges = config->num_hyperedges;
    hg->total_elements = config->num_nodes + config->num_hyperedges;

    hg->node_degrees = (int*)calloc(hg->total_elements, sizeof(int));
    hg->hyperedge_sizes = (int*)calloc(hg->num_hyperedges, sizeof(int));
    hg->visited = (int*)calloc(hg->total_elements, sizeof(int));
    hg->damaged = (int*)calloc(hg->total_elements, sizeof(int));
    hg->temp_state = (int*)calloc(hg->total_elements, sizeof(int));
    hg->cluster_sizes = (int*)calloc(hg->total_elements, sizeof(int));

    if (!hg->node_degrees || !hg->hyperedge_sizes || !hg->visited ||
        !hg->damaged || !hg->temp_state || !hg->cluster_sizes) {
        fprintf(stderr, "Error: Memory allocation failed for hypergraph arrays\n");
        free_hypergraph(hg);
        return NULL;
    }

    hg->adjacency_lists = (int**)calloc(hg->total_elements, sizeof(int*));
    hg->weights = (float**)calloc(hg->total_elements, sizeof(float*));

    if (!hg->adjacency_lists || !hg->weights) {
        fprintf(stderr, "Error: Memory allocation failed for adjacency lists\n");
        free_hypergraph(hg);
        return NULL;
    }


    return hg;
}

void free_hypergraph(hypergraph_t* hg) {
    if (!hg) return;

    // Free adjacency lists and weights
    if (hg->adjacency_lists) {
        for (int i = 0; i < hg->total_elements; i++) {
            if (hg->adjacency_lists[i]) {
                free(hg->adjacency_lists[i]);
            }
        }
        free(hg->adjacency_lists);
    }

    if (hg->weights) {
        for (int i = 0; i < hg->total_elements; i++) {
            if (hg->weights[i]) {
                free(hg->weights[i]);
            }
        }
        free(hg->weights);
    }

    free(hg->node_degrees);
    free(hg->hyperedge_sizes);
    free(hg->visited);
    free(hg->damaged);
    free(hg->temp_state);
    free(hg->cluster_sizes);

    // Free the structure itself
    free(hg);
}

bool generate_random_hypergraph(hypergraph_t* hg, const config_t* config) {
    if (!hg || !config) return false;

    int N = config->num_nodes;
    int Q = config->num_hyperedges;
    float mav = config->avg_membership;

    // First, determine the size (number of nodes) for each hyperedge
    for (int mu = 0; mu < Q; mu++) {
        // Start with base size of 2
        hg->hyperedge_sizes[mu] = 2;

        // Probabilistically add more nodes based on mav
        for (int i = 0; i < N; i++) {
            if (random_double() < (mav - 2.0) / N) {
                hg->hyperedge_sizes[mu]++;
            }
        }
    }

    for (int i = 0; i < N + Q; i++) {
         (i < N), 
         (i >= N), 
        int max_connections = (i < N) ? N : hg->hyperedge_sizes[i - N];

        hg->adjacency_lists[i] = (int*)calloc(max_connections, sizeof(int));
        hg->weights[i] = (float*)calloc(max_connections, sizeof(float));

        if (!hg->adjacency_lists[i] || !hg->weights[i]) {
            fprintf(stderr, "Error: Memory allocation failed for adjacency list %d\n", i);
            return false;
        }
    }

    for (int mu = 0; mu < Q; mu++) {
        int hyperedge_id = N + mu;
        int* included_nodes = (int*)calloc(hg->hyperedge_sizes[mu], sizeof(int));

        if (!included_nodes) {
            fprintf(stderr, "Error: Memory allocation failed for temporary node list\n");
            return false;
        }

        for (int n = 0; n < hg->hyperedge_sizes[mu]; n++) {
            int node_id;
            bool already_included;

            do {
                node_id = (int)((float)N * random_double());

                already_included = false;
                for (int nn = 0; nn < n; nn++) {
                    if (included_nodes[nn] == node_id) {
                        already_included = true;
                        break;
                    }
                }
            } while (already_included);

            included_nodes[n] = node_id;

            hg->adjacency_lists[node_id][hg->node_degrees[node_id]] = hyperedge_id;
            hg->adjacency_lists[hyperedge_id][n] = node_id;

            hg->weights[node_id][hg->node_degrees[node_id]] = (int)(2 * random_double());
            hg->weights[hyperedge_id][n] = (int)(2 * random_double());

            hg->node_degrees[node_id]++;
        }

        hg->node_degrees[hyperedge_id] = hg->hyperedge_sizes[mu];

        free(included_nodes);
    }

    for (int i = 0; i < N; i++) {
        int degree = hg->node_degrees[i];

        if (degree < N) {
            int* new_adjacency = (int*)realloc(hg->adjacency_lists[i], degree * sizeof(int));
            float* new_weights = (float*)realloc(hg->weights[i], degree * sizeof(float));

            if (degree > 0 && (!new_adjacency || !new_weights)) {
                fprintf(stderr, "Error: Memory reallocation failed for node %d\n", i);
                return false;
            }

            if (new_adjacency) hg->adjacency_lists[i] = new_adjacency;
            if (new_weights) hg->weights[i] = new_weights;
        }
    }

    return true;
}

int find_connected_component(hypergraph_t* hg, int start_node, int cluster_id, const config_t* config) {
    int cluster_size = 0;
    int s_core = config->s_core;

    hg->visited[start_node] = cluster_id;

    if (start_node < hg->num_nodes) {
        cluster_size++;
    }

    for (int n = 0; n < hg->node_degrees[start_node]; n++) {
        int neighbor = hg->adjacency_lists[start_node][n];

        if (neighbor < hg->num_nodes) {
            if (hg->visited[neighbor] == 0 && hg->damaged[neighbor] == 1) {
                int sub_cluster_size = find_connected_component(hg, neighbor, cluster_id, config);
                cluster_size += sub_cluster_size;
            }
        } else {
            if (hg->node_degrees[neighbor] >= s_core) {
                int all_undamaged = 1;

                for (int j = 0; j < hg->node_degrees[neighbor]; j++) {
                    int node = hg->adjacency_lists[neighbor][j];
                    all_undamaged &= hg->damaged[node];
                }

                if (all_undamaged && hg->visited[neighbor] == 0) {
                    int sub_cluster_size = find_connected_component(hg, neighbor, cluster_id, config);
                    cluster_size += sub_cluster_size;
                }
            }
        }
    }

    return cluster_size;
}


bool run_percolation_simulation(hypergraph_t* hg, const config_t* config, float* results) {
    int N = hg->num_nodes;
    int num_points = 100;  // Number of data points (pruning percentages)
    float* damage_thresholds = (float*)calloc(N, sizeof(float));

    if (!damage_thresholds) {
        fprintf(stderr, "Error: Memory allocation failed for damage thresholds\n");
        return false;
    }

    for (int nc = 0; nc < num_points; nc++) {
        results[nc] = 0.0;
    }

    for (int nrun = 0; nrun < config->num_simulations; nrun++) {
        if (config->verbose && nrun % 10 == 0) {
            print_progress(nrun, config->num_simulations);
        }

        for (int i = 0; i < N; i++) {
            damage_thresholds[i] = random_double();
        }

        for (int nc = 0; nc < num_points; nc++) {
            float pruning_threshold = nc * 0.01;  // 0.00, 0.01, ..., 0.99

            for (int i = 0; i < N; i++) {
                hg->damaged[i] = (damage_thresholds[i] >= pruning_threshold) ? 1 : 0;
            }

            bool changed;
            do {
                changed = false;

                for (int i = 0; i < N; i++) {
                    hg->temp_state[i] = hg->damaged[i];
                }

                for (int i = 0; i < N; i++) {
                    if (hg->damaged[i] == 1) {
                        int valid_connections = 0;

                        for (int n = 0; n < hg->node_degrees[i]; n++) {
                            int hyperedge = hg->adjacency_lists[i][n];

                            if (hg->node_degrees[hyperedge] >= config->s_core) {
                                int all_valid = 1;

                                for (int j = 0; j < hg->node_degrees[hyperedge]; j++) {
                                    int node = hg->adjacency_lists[hyperedge][j];
                                    all_valid &= hg->damaged[node];
                                }

                                if (all_valid) {
                                    valid_connections++;
                                }
                            }
                        }

                        if (valid_connections < config->k_core) {
                            hg->temp_state[i] = 0;
                            changed = true;
                        }
                    }
                }

                for (int i = 0; i < N; i++) {
                    hg->damaged[i] = hg->temp_state[i];
                }

            } while (changed);

            for (int i = 0; i < hg->total_elements; i++) {
                hg->visited[i] = 0;
                hg->cluster_sizes[i] = 0;
            }

            int num_clusters = 0;
            int largest_cluster_size = 0;

            for (int i = 0; i < N; i++) {
                if (hg->visited[i] == 0 && hg->damaged[i] == 1) {
                    num_clusters++;
                    int cluster_size = find_connected_component(hg, i, num_clusters, config);
                    hg->cluster_sizes[num_clusters] = cluster_size;

                    if (cluster_size > largest_cluster_size) {
                        largest_cluster_size = cluster_size;
                        // We don't need to keep track of the cluster ID, just its size
                    }
                }
            }

            results[nc] += (float)largest_cluster_size / (float)N;
        }
    }

    for (int nc = 0; nc < num_points; nc++) {
        results[nc] /= config->num_simulations;
    }

    if (config->verbose) {
        print_progress(config->num_simulations, config->num_simulations);
        printf("\n");
    }

    free(damage_thresholds);
    return true;
}

bool write_results_to_file(const float* results, int num_points, const config_t* config) {
    FILE* fp = fopen(config->output_filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: Could not open file '%s' for writing: %s\n",
                config->output_filename, strerror(errno));
        return false;
    }

    fprintf(fp, "# Hypergraph K-Core Percolation Simulation Results\n");
    fprintf(fp, "# Parameters: N=%d, Q=%d, k=%d, s=%d, m_avg=%.2f, runs=%d\n",
            config->num_nodes, config->num_hyperedges, config->k_core,
            config->s_core, config->avg_membership, config->num_simulations);
    fprintf(fp, "# Format: 1-p, GCS/N\n");
    fprintf(fp, "# Where p is the pruning probability and GCS is the Giant Component Size\n");

    for (int nc = 0; nc < num_points; nc++) {
        float pruning_probability = nc * 0.01;
        fprintf(fp, "%.6f %.6f\n", 1.0 - pruning_probability, results[nc]);
    }

    fclose(fp);
    return true;
}

void print_progress(int current, int total) {
    const int bar_width = 50;
    float progress = (float)current / total;
    int pos = bar_width * progress;

    printf("[");
    for (int i = 0; i < bar_width; i++) {
        if (i < pos) printf("=");
        else if (i == pos) printf(">");
        else printf(" ");
    }
    printf("] %d%%\r", (int)(progress * 100.0));
    fflush(stdout);
}

int main(int argc, char** argv) {
    config_t config;
    hypergraph_t* hypergraph = NULL;
    float* results = NULL;
    bool success = true;

    if (!parse_arguments(argc, argv, &config)) {
        return EXIT_FAILURE;
    }

    if (config.random_seed == 0) {
        init_random((unsigned long)time(NULL));
    } else {
        init_random(config.random_seed);
    }

    if (config.verbose) {
        printf("Hypergraph K-Core Percolation Simulator\n");
        printf("=======================================\n");
        printf("Nodes: %d\n", config.num_nodes);
        printf("Hyperedges: %d\n", config.num_hyperedges);
        printf("k-core parameter: %d\n", config.k_core);
        printf("s-core parameter: %d\n", config.s_core);
        printf("Average membership: %.2f\n", config.avg_membership);
        printf("Monte Carlo simulations: %d\n", config.num_simulations);
        printf("Output file: %s\n", config.output_filename);
        printf("Random seed: %lu\n", (unsigned long)(config.random_seed == 0 ? time(NULL) : config.random_seed));
        printf("=======================================\n\n");
    }

    // Allocate memory for results
    results = (float*)calloc(100, sizeof(float));
    if (!results) {
        fprintf(stderr, "Error: Memory allocation failed for results array\n");
        return EXIT_FAILURE;
    }

    // Create hypergraph
    if (config.verbose) printf("Creating hypergraph structure...\n");
    hypergraph = create_hypergraph(&config);
    if (!hypergraph) {
        free(results);
        return EXIT_FAILURE;
    }

    // Generate random hypergraph
    if (config.verbose) printf("Generating random hypergraph...\n");
    if (!generate_random_hypergraph(hypergraph, &config)) {
        fprintf(stderr, "Error: Failed to generate random hypergraph\n");
        free_hypergraph(hypergraph);
        free(results);
        return EXIT_FAILURE;
    }

    // Run percolation simulation
    if (config.verbose) printf("Running percolation simulation...\n");
    if (!run_percolation_simulation(hypergraph, &config, results)) {
        fprintf(stderr, "Error: Failed to run percolation simulation\n");
        free_hypergraph(hypergraph);
        free(results);
        return EXIT_FAILURE;
    }

    // Write results to file
    if (config.verbose) printf("Writing results to file...\n");
    if (!write_results_to_file(results, 100, &config)) {
        fprintf(stderr, "Error: Failed to write results to file\n");
        success = false;
    }

    // Clean up
    free_hypergraph(hypergraph);
    free(results);

    if (config.verbose) {
        if (success) {
            printf("Simulation completed successfully!\n");
            printf("Results written to '%s'\n", config.output_filename);
        } else {
            printf("Simulation completed with errors.\n");
        }
    }

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
