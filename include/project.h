#ifndef PROJECT_H
#define PROJECT_H

#include <stddef.h>

#define MAX_SEQ_LEN 1024
#define MAX_SEQ_COUNT 100

typedef struct {
    char id[50];
    char sequence[MAX_SEQ_LEN];
} Sequence;

typedef struct {
    size_t n; // number of sequences
    Sequence sequences[MAX_SEQ_COUNT];
} SequenceCollection;

void perform_fft(const char *sequence, double *out_fft, size_t length);
void pairwise_alignment(const char *seq1, const char *seq2, double *score);
void compute_distance_matrix(const SequenceCollection *seq_collection, double **distance_matrix);

#endif
