#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "project.h"

/*
 * perform_fft:
 * Dummy implementation of FFT transformation.
 * convert each character to its double value.
 */
void perform_fft(const char *sequence, double *out_fft, size_t length) {
    for (size_t i = 0; i < length; i++) {
        out_fft[i] = (double)sequence[i];
    }
}

/*
 * pairwise_alignment:
 * A dummy implementation of pairwise alignment that counts mismatches.
 */
void pairwise_alignment(const char *seq1, const char *seq2, double *score) {
    size_t len1 = strlen(seq1);
    size_t len2 = strlen(seq2);
    size_t min_len = (len1 < len2) ? len1 : len2;
    double mismatches = 0.0;
    for (size_t i = 0; i < min_len; i++) {
        if (seq1[i] != seq2[i]) {
            mismatches += 1.0;
        }
    }
    // Add a penalty for length difference.
    mismatches += fabs((double)len1 - (double)len2);
    *score = mismatches;
}

/*
 * compute_distance_matrix:
 * For each pair of sequences, perform a dummy alignment and fill in the distance matrix.
 */
void compute_distance_matrix(const SequenceCollection *seq_collection, double **distance_matrix) {
    size_t n = seq_collection->n;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            if (i == j) {
                distance_matrix[i][j] = 0.0;
            } else {
                double score;
                pairwise_alignment(seq_collection->sequences[i].sequence, 
                                   seq_collection->sequences[j].sequence, 
                                   &score);
                // In this demo, the alignment score directly represents the distance.
                distance_matrix[i][j] = score;
            }
        }
    }
}
