#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include "utils.h"

SequenceCollection* read_fasta(const char* filename) {
    FILE* file = fopen(filename, "r");
    if(!file) {
        perror("Failed to open file");
        return NULL;
    }

    SequenceCollection* collection = malloc(sizeof(SequenceCollection));
    collection->count = 0;
    collection->sequences = NULL;
    collection->fft_length = 0;

    char buffer[1024];
    char* current_seq = NULL;
    size_t seq_len = 0;

    while(fgets(buffer, sizeof(buffer), file)) {
        if(buffer[0] == '>') {
            // Finalize previous sequence
            if(current_seq) {
                collection->sequences = realloc(collection->sequences, 
                    (collection->count+1) * sizeof(Sequence));
                collection->sequences[collection->count].sequence = current_seq;
                size_t encoded_len;
                collection->sequences[collection->count].encoded = 
                    encode_sequence(current_seq, &encoded_len);
                collection->sequences[collection->count].length = encoded_len;
                collection->count++;
            }
            current_seq = malloc(1024);
            current_seq[0] = '\0';
            seq_len = 0;
        } else {
            // Remove newline characters
            char* nl = strchr(buffer, '\n');
            if(nl) *nl = '\0';
            size_t len = strlen(buffer);
            current_seq = realloc(current_seq, seq_len + len + 1);
            strcpy(current_seq + seq_len, buffer);
            seq_len += len;
        }
    }

    if(current_seq) {
        collection->sequences = realloc(collection->sequences, 
            (collection->count+1) * sizeof(Sequence));
        collection->sequences[collection->count].sequence = current_seq;
        size_t encoded_len;
        collection->sequences[collection->count].encoded = 
            encode_sequence(current_seq, &encoded_len);
        collection->sequences[collection->count].length = encoded_len;
        collection->count++;
    }

    fclose(file);
    return collection;
}

double** create_matrix(size_t n) {
    double** matrix = malloc(n * sizeof(double*));
    for(size_t i=0; i<n; i++) {
        matrix[i] = malloc(n * sizeof(double));
        memset(matrix[i], 0, n * sizeof(double));
    }
    return matrix;
}

void start_timer(Timer* t, const char* name) {
    t->name = name;
    t->start = clock();
    printf("[START] %s\n", name);
}

void stop_timer(Timer* t) {
    t->end = clock();
    printf("[END] %s\n", t->name);
}

void print_elapsed(Timer t) {
    double elapsed = (double)(t.end - t.start) * 1000.0 / CLOCKS_PER_SEC;
    printf("[TIMER] %s: %.2f ms\n\n", t.name, elapsed);
}

void print_distance_matrix(double** matrix, size_t n) {
    printf("\nDistance Matrix:\n");
    for(size_t i=0; i<n; i++) {
        for(size_t j=0; j<n; j++) {
            printf("%6.3f\t", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_encoded_sequence(const double* encoded, size_t len) {
    printf("\nEncoded Sequence:\n");
    for(size_t i=0; i<len; i++) {
        printf("Position %zu: ", i);
        for(size_t k=0; k<ENCODING_DIM; k++) {
            printf("%6.2f ", encoded[i*ENCODING_DIM + k]);
        }
        printf("\n");
    }
}

void print_fft_data(const complex_t** fft_data, size_t n_padded, size_t dim) {
    printf("\nFFT Data (Dimension %zu):\n", dim);
    for(size_t i=0; i<n_padded; i++) {
        printf("%6.2f+%6.2fi ", fft_data[dim][i].real, fft_data[dim][i].imag);
        if((i+1) % 4 == 0) printf("\n");
    }
    printf("\n");
}

void print_tree(const UPGMA_Node* node, int level) {
    if(!node) return;
    
    for(int i=0; i<level; i++) printf("  ");
    
    if(node->id >= 0) {
        printf("Leaf %d (Height: %.3f)\n", node->id, node->height);
    } else {
        printf("Node (Size: %zu, Height: %.3f)\n", node->size, node->height);
        print_tree(node->left, level+1);
        print_tree(node->right, level+1);
    }
}

void print_guide_tree(const UPGMA_Node* root) {
    printf("\nGuide Tree Structure:\n");
    print_tree(root, 0);
}

void free_matrix(double** matrix, size_t n) {
    for(size_t i=0; i<n; i++) free(matrix[i]);
    free(matrix);
}

void free_tree(UPGMA_Node* node) {
    if(!node) return;
    free_tree(node->left);
    free_tree(node->right);
    free(node->alignment);
    free(node);
}

void free_collection(SequenceCollection* collection) {
    for(size_t i=0; i<collection->count; i++) {
        free(collection->sequences[i].sequence);
        free(collection->sequences[i].encoded);
        for(size_t k=0; k<ENCODING_DIM; k++)
            free(collection->sequences[i].fft_data[k]);
        free(collection->sequences[i].fft_data);
    }
    free(collection->sequences);
    free(collection);
}