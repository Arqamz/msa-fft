#ifndef PROJECT_H
#define PROJECT_H

#include <stdio.h>
#include <stddef.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>

#define ENCODING_DIM 5

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

typedef struct {
    double real;
    double imag;
} complex_t;

typedef struct {
    char* sequence;
    double* encoded;
    size_t length;
    complex_t** fft_data;
} Sequence;

typedef struct {
    Sequence* sequences;
    size_t count;
    size_t fft_length;
} SequenceCollection;

typedef struct UPGMA_Node {
    int id;
    struct UPGMA_Node* left;
    struct UPGMA_Node* right;
    double height;
    size_t size;
    char* alignment;
} UPGMA_Node;

/* Encoding functions */
double* encode_sequence(const char* sequence, size_t* out_length);

/* FFT functions */
size_t next_pow2(size_t n);
void fft(complex_t* data, size_t n);
void inverse_fft(complex_t* data, size_t n);
void perform_fft(const double* encoded_seq, size_t len, size_t n_padded, complex_t*** fft_result);
void perform_fft_analysis(SequenceCollection* collection);
double compute_max_cross_correlation(complex_t** fft_a, complex_t** fft_b, size_t n_padded);

/* Distance matrix */
double** compute_fft_distance_matrix(SequenceCollection* collection);

/* Guide tree construction */
UPGMA_Node* build_guide_tree(double** distance_matrix, size_t n);

/* Alignment functions */
char* needleman_wunsch_align(const char* seq1, const char* seq2);
void progressive_alignment(SequenceCollection* collection, UPGMA_Node* tree);

/* Utility functions */
double** create_matrix(size_t n);
void free_matrix(double** matrix, size_t n);
SequenceCollection* read_fasta(const char* filename);
void free_tree(UPGMA_Node* node);
void free_collection(SequenceCollection* collection);

#endif /* PROJECT_H */