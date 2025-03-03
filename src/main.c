#include <stdio.h>
#include <stdlib.h>
#include "project.h"
#include "utils.h"

int main() {
    // Load sequences (currently hard-coded)
    SequenceCollection seq_collection;
    load_sequences(&seq_collection);

    size_t n = seq_collection.n;
    
    // Allocate memory for the distance matrix (n x n)
    double **distance_matrix = malloc(n * sizeof(double *));
    for (size_t i = 0; i < n; i++) {
        distance_matrix[i] = malloc(n * sizeof(double));
    }

    // Compute the distance matrix based on pairwise alignment scores.
    compute_distance_matrix(&seq_collection, distance_matrix);
    print_distance_matrix(distance_matrix, n);

    // Free allocated memory.
    for (size_t i = 0; i < n; i++) {
        free(distance_matrix[i]);
    }
    free(distance_matrix);

    return 0;
}
