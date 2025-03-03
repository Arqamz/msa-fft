#ifndef UTILS_H
#define UTILS_H

#include "project.h"
#include <time.h>

typedef struct {
    clock_t start;
    clock_t end;
    const char* name;
} Timer;

// Timing functions
void start_timer(Timer* t, const char* name);
void stop_timer(Timer* t);
void print_elapsed(Timer t);

// Matrix operations
void print_distance_matrix(double** matrix, size_t n);
void print_encoded_sequence(const double* encoded, size_t len);
void print_fft_data(const complex_t** fft_data, size_t n_padded, size_t dim);

// Tree debugging
void print_tree(const UPGMA_Node* node, int level);
void print_guide_tree(const UPGMA_Node* root);

// Existing functions
SequenceCollection* read_fasta(const char* filename);
double** create_matrix(size_t n);
void free_matrix(double** matrix, size_t n);
void free_tree(UPGMA_Node* node);
void free_collection(SequenceCollection* collection);

#endif