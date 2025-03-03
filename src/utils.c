#include <stdio.h>
#include <string.h>
#include "utils.h"

/*
 * print_distance_matrix:
 * Prints the distance matrix in a formatted way.
 */
void print_distance_matrix(double **matrix, size_t n) {
    printf("Distance Matrix:\n");
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            printf("%6.2f ", matrix[i][j]);
        }
        printf("\n");
    }
}

/*
 * load_sequences:
 * Loads a set of hard-coded sequences.
 * For the full implementation, this function will read from a file.
 */
void load_sequences(SequenceCollection *collection) {
    collection->n = 4;
    
    strcpy(collection->sequences[0].id, "Seq1");
    strcpy(collection->sequences[0].sequence, "MKTLLILTCLVAVALARPKAQQL");
    
    strcpy(collection->sequences[1].id, "Seq2");
    strcpy(collection->sequences[1].sequence, "MKTVLILTCLVALAKPKAQQL");
    
    strcpy(collection->sequences[2].id, "Seq3");
    strcpy(collection->sequences[2].sequence, "MKTLLILACLVALARKAQQL");
    
    strcpy(collection->sequences[3].id, "Seq4");
    strcpy(collection->sequences[3].sequence, "MKTLLILTCLVALAKPQQL");
}
