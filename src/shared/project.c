#include "project.h"
#include "utils.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/*
 * encode_sequence:
 * This function converts an amino acid sequence into a numerical representation.
 * Each amino acid is mapped to a 5-dimensional vector of physicochemical properties:
 *   [hydrophobicity, volume, polarity, isoelectric point, molecular weight]
 *
 * After mapping, each column (property) is normalized to have zero mean and unit variance.
 *
 * The function returns a pointer to a dynamically allocated array of doubles with size:
 *     strlen(sequence) * ENCODING_DIM.
 * The number of residues (sequence length) is returned via the out_length pointer.
 */
double* encode_sequence(const char* sequence, size_t *out_length) {
    size_t len = strlen(sequence);
    *out_length = len;
    double *encoded = malloc(len * ENCODING_DIM * sizeof(double));
    if (!encoded) {
        perror("Failed to allocate memory for encoded sequence");
        exit(EXIT_FAILURE);
    }
    
    // Define a mapping for each amino acid to a 5-dimensional vector.
    typedef struct {
        char aa;
        double vec[ENCODING_DIM];
    } AAEncoding;
    
    // The following values are illustrative and based on selected physicochemical properties.
    static const AAEncoding encoding_table[] = {
        { 'A', { 1.8,   88.6,   8.1,  6.00,  89.1 } },
        { 'R', { -4.5, 173.4,  10.5, 10.76, 174.2 } },
        { 'N', { -3.5, 114.1,  11.6,  5.41, 132.1 } },
        { 'D', { -3.5, 111.1,  13.0,  2.77, 133.1 } },
        { 'C', { 2.5,  108.5,   5.5,  5.07, 121.2 } },
        { 'E', { -3.5, 138.4,  12.3,  3.22, 147.1 } },
        { 'Q', { -3.5, 143.8,  10.5,  5.65, 146.1 } },
        { 'G', { -0.4,  60.1,   9.0,  5.97,  75.1 } },
        { 'H', { -3.2, 153.2,  10.4,  7.59, 155.2 } },
        { 'I', { 4.5,  166.7,   5.2,  6.02, 131.2 } },
        { 'L', { 3.8,  166.7,   4.9,  6.01, 131.2 } },
        { 'K', { -3.9, 168.6,  11.3,  9.74, 146.2 } },
        { 'M', { 1.9,  162.9,   5.7,  5.74, 149.2 } },
        { 'F', { 2.8,  189.9,   5.0,  5.48, 165.2 } },
        { 'P', { -1.6, 112.7,   8.0,  6.30, 115.1 } },
        { 'S', { -0.8,  89.0,   9.2,  5.68, 105.1 } },
        { 'T', { -0.7, 116.1,   8.6,  5.60, 119.1 } },
        { 'W', { -0.9, 227.8,   5.4,  5.89, 204.2 } },
        { 'Y', { -1.3, 193.6,   6.2,  5.66, 181.2 } },
        { 'V', { 4.2,  140.0,   5.9,  6.00, 117.1 } }
    };
    const size_t table_size = sizeof(encoding_table) / sizeof(AAEncoding);
    
    // Map each amino acid in the input sequence to its corresponding vector.
    for (size_t i = 0; i < len; i++) {
        char aa = toupper(sequence[i]);
        int found = 0;
        for (size_t j = 0; j < table_size; j++) {
            if (encoding_table[j].aa == aa) {
                for (size_t k = 0; k < ENCODING_DIM; k++) {
                    encoded[i * ENCODING_DIM + k] = encoding_table[j].vec[k];
                }
                found = 1;
                break;
            }
        }
        // If the amino acid is not in the table, assign a zero vector.
        if (!found) {
            for (size_t k = 0; k < ENCODING_DIM; k++) {
                encoded[i * ENCODING_DIM + k] = 0.0;
            }
        }
    }
    
    // Normalize each property (each column) across the sequence.
    for (size_t k = 0; k < ENCODING_DIM; k++) {
        double sum = 0.0;
        for (size_t i = 0; i < len; i++) {
            sum += encoded[i * ENCODING_DIM + k];
        }
        double mean = sum / len;
        
        double variance = 0.0;
        for (size_t i = 0; i < len; i++) {
            double diff = encoded[i * ENCODING_DIM + k] - mean;
            variance += diff * diff;
        }
        double std = (len > 1) ? sqrt(variance / (len - 1)) : 1.0;
        if (std == 0.0) std = 1.0;
        
        // Normalize the k-th column.
        for (size_t i = 0; i < len; i++) {
            encoded[i * ENCODING_DIM + k] = (encoded[i * ENCODING_DIM + k] - mean) / std;
        }
    }
    
    return encoded;
}

size_t next_pow2(size_t n) {
    if (n == 0) return 1;
    size_t p = 1;
    while (p < n) p <<= 1;
    return p;
}

void fft(complex_t *data, size_t n) {
    // Bit-reversal permutation
    for (size_t i = 0, j = 0; i < n; i++) {
        if (j > i) {
            complex_t temp = data[i];
            data[i] = data[j];
            data[j] = temp;
        }
        size_t m = n >> 1;
        while (m >= 1 && j >= m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }

    // Iterative FFT
    for (size_t s = 1; s <= (size_t)log2(n); s++) {
        size_t m = 1 << s;
        size_t mhalf = m >> 1;
        for (size_t k = 0; k < n; k += m) {
            for (size_t j = 0; j < mhalf; j++) {
                double angle = -2.0 * M_PI * j / m;
                complex_t w = { cos(angle), sin(angle) };
                complex_t t = {
                    w.real * data[k + j + mhalf].real - w.imag * data[k + j + mhalf].imag,
                    w.real * data[k + j + mhalf].imag + w.imag * data[k + j + mhalf].real
                };
                complex_t even = data[k + j];
                data[k + j].real = even.real + t.real;
                data[k + j].imag = even.imag + t.imag;
                data[k + j + mhalf].real = even.real - t.real;
                data[k + j + mhalf].imag = even.imag - t.imag;
            }
        }
    }
}

void inverse_fft(complex_t *data, size_t n) {
    // Conjugate the data
    for (size_t i = 0; i < n; i++) data[i].imag = -data[i].imag;
    
    // Compute FFT
    fft(data, n);
    
    // Conjugate again and scale
    double scale = 1.0 / n;
    for (size_t i = 0; i < n; i++) {
        data[i].real *= scale;
        data[i].imag = -data[i].imag * scale;
    }
}

// FFT processing for sequences
void perform_fft_analysis(SequenceCollection* collection) {
    size_t max_len = 0;
    for(size_t i=0; i<collection->count; i++) {
        if(collection->sequences[i].length > max_len)
            max_len = collection->sequences[i].length;
    }
    
    size_t n_padded = next_pow2(2 * max_len - 1);
    collection->fft_length = n_padded;
    
    for(size_t i=0; i<collection->count; i++) {
        complex_t*** fft_result = &collection->sequences[i].fft_data;
        perform_fft(collection->sequences[i].encoded, 
                   collection->sequences[i].length,
                   n_padded, fft_result);
    }
}

void perform_fft(const double* encoded_seq, size_t len, size_t n_padded, complex_t ***fft_result) {
    *fft_result = malloc(ENCODING_DIM * sizeof(complex_t *));
    for(size_t k=0; k<ENCODING_DIM; k++) {
        (*fft_result)[k] = malloc(n_padded * sizeof(complex_t));
        for(size_t i=0; i<n_padded; i++) {
            if(i < len) {
                (*fft_result)[k][i].real = encoded_seq[i * ENCODING_DIM + k];
            } else {
                (*fft_result)[k][i].real = 0.0;
            }
            (*fft_result)[k][i].imag = 0.0;
        }
        fft((*fft_result)[k], n_padded);
    }
}

// Distance matrix calculation
double compute_max_cross_correlation(complex_t **fft_a, complex_t **fft_b, size_t n_padded) {
    complex_t *total_cross = calloc(n_padded, sizeof(complex_t));
    
    for(size_t k=0; k<ENCODING_DIM; k++) {
        complex_t *cross = malloc(n_padded * sizeof(complex_t));
        for(size_t m=0; m<n_padded; m++) {
            complex_t b_conj = {fft_b[k][m].real, -fft_b[k][m].imag};
            cross[m].real = fft_a[k][m].real * b_conj.real - fft_a[k][m].imag * b_conj.imag;
            cross[m].imag = fft_a[k][m].real * b_conj.imag + fft_a[k][m].imag * b_conj.real;
        }
        inverse_fft(cross, n_padded);
        
        for(size_t m=0; m<n_padded; m++) {
            total_cross[m].real += cross[m].real;
        }
        free(cross);
    }

    double max_corr = -DBL_MAX;
    for(size_t m=0; m<n_padded; m++) {
        if(total_cross[m].real > max_corr) max_corr = total_cross[m].real;
    }
    free(total_cross);
    return max_corr;
}

double** compute_fft_distance_matrix(SequenceCollection* collection) {
    perform_fft_analysis(collection);
    double** matrix = create_matrix(collection->count);
    
    for(size_t i=0; i<collection->count; i++) {
        for(size_t j=i+1; j<collection->count; j++) {
            double corr = compute_max_cross_correlation(
                collection->sequences[i].fft_data,
                collection->sequences[j].fft_data,
                collection->fft_length
            );
            matrix[i][j] = matrix[j][i] = 1.0 - (corr + 1.0)/2.0;
        }
    }
    return matrix;
}

// UPGMA Tree Construction
UPGMA_Node* build_guide_tree(double** distance_matrix, size_t n) {
    UPGMA_Node** nodes = malloc(n * sizeof(UPGMA_Node*));
    size_t* sizes = malloc(n * sizeof(size_t));
    
    // Initialize leaves
    for(size_t i = 0; i < n; i++) {
        nodes[i] = malloc(sizeof(UPGMA_Node));
        nodes[i]->id = i;
        nodes[i]->left = nodes[i]->right = NULL;
        nodes[i]->height = 0.0;
        nodes[i]->size = 1;
        sizes[i] = 1;
    }

    size_t remaining = n;
    size_t a = 0, b = 0;

    // Initialize a and b to valid indices to avoid uninitialized usage warning
    int found = 0; // A flag to indicate if a valid pair (a, b) is found
    
    while(remaining > 1) {
        // Find closest clusters
        double min_dist = INFINITY;
        
        for(size_t i = 0; i < n; i++) {
            if(!nodes[i]) continue;
            for(size_t j = i + 1; j < n; j++) {
                if(!nodes[j]) continue;
                if(distance_matrix[i][j] < min_dist) {
                    min_dist = distance_matrix[i][j];
                    a = i;
                    b = j;
                    found = 1; // Mark that we found a valid pair
                }
            }
        }

        if (!found) {
            // If no valid pair was found, break the loop
            break;
        }

        // Create new cluster
        UPGMA_Node* new_node = malloc(sizeof(UPGMA_Node));
        new_node->id = -1;
        new_node->left = nodes[a];
        new_node->right = nodes[b];
        new_node->height = min_dist / 2.0;
        new_node->size = sizes[a] + sizes[b];

        // Update distance matrix
        for(size_t i = 0; i < n; i++) {
            if(i == a || i == b || !nodes[i]) continue;
            distance_matrix[a][i] = distance_matrix[i][a] = 
                (sizes[a] * distance_matrix[a][i] + sizes[b] * distance_matrix[b][i]) / 
                (sizes[a] + sizes[b]);
        }

        sizes[a] += sizes[b];
        nodes[b] = NULL;
        nodes[a] = new_node;
        remaining--;
    }

    UPGMA_Node* root = NULL;
    if (a < n && nodes[a]) {
        root = nodes[a];
    }

    free(nodes);
    free(sizes);
    return root;
}

// Progressive Alignment Core
char* needleman_wunsch_align(const char* seq1, const char* seq2) {
    const int MATCH = 1;
    const int MISMATCH = -1;
    const int GAP_PENALTY = -2;

    size_t len1 = strlen(seq1);
    size_t len2 = strlen(seq2);

    // Create DP table
    int** score = malloc((len1+1) * sizeof(int*));
    char** trace = malloc((len1+1) * sizeof(char*));
    
    for(size_t i=0; i<=len1; i++) {
        score[i] = malloc((len2+1) * sizeof(int));
        trace[i] = malloc((len2+1) * sizeof(char));
    }

    // Initialize matrices
    for(size_t i=0; i<=len1; i++) {
        score[i][0] = i * GAP_PENALTY;
        trace[i][0] = 'U';
    }
    for(size_t j=0; j<=len2; j++) {
        score[0][j] = j * GAP_PENALTY;
        trace[0][j] = 'L';
    }

    // Fill matrices
    for(size_t i=1; i<=len1; i++) {
        for(size_t j=1; j<=len2; j++) {
            int match_score = (toupper(seq1[i-1]) == toupper(seq2[j-1])) ? MATCH : MISMATCH;
            int diag = score[i-1][j-1] + match_score;
            int up = score[i-1][j] + GAP_PENALTY;
            int left = score[i][j-1] + GAP_PENALTY;

            if(diag >= up && diag >= left) {
                score[i][j] = diag;
                trace[i][j] = 'D';
            } else if(up >= left) {
                score[i][j] = up;
                trace[i][j] = 'U';
            } else {
                score[i][j] = left;
                trace[i][j] = 'L';
            }
        }
    }

    // Traceback
    size_t i = len1, j = len2;
    size_t len = len1 + len2;
    char* aligned1 = malloc(len+1);
    char* aligned2 = malloc(len+1);
    size_t pos = 0;

    while(i > 0 || j > 0) {
        if(trace[i][j] == 'D') {
            aligned1[pos] = seq1[i-1];
            aligned2[pos] = seq2[j-1];
            i--;
            j--;
        } else if(trace[i][j] == 'U') {
            aligned1[pos] = seq1[i-1];
            aligned2[pos] = '-';
            i--;
        } else {
            aligned1[pos] = '-';
            aligned2[pos] = seq2[j-1];
            j--;
        }
        pos++;
    }
    aligned1[pos] = '\0';
    aligned2[pos] = '\0';

    // Reverse the strings
    size_t start = 0, end = pos-1;
    while(start < end) {
        char tmp = aligned1[start];
        aligned1[start] = aligned1[end];
        aligned1[end] = tmp;
        tmp = aligned2[start];
        aligned2[start] = aligned2[end];
        aligned2[end] = tmp;
        start++;
        end--;
    }

    // Merge aligned sequences
    char* merged = malloc(2 * pos + 1);
    snprintf(merged, 2 * pos + 1, "%s\n%s", aligned1, aligned2);

    // Cleanup
    free(aligned1);
    free(aligned2);
    for(size_t i=0; i<=len1; i++) {
        free(score[i]);
        free(trace[i]);
    }
    free(score);
    free(trace);

    return merged;
}

// Progressive alignment driver
void progressive_alignment(SequenceCollection* collection, UPGMA_Node* tree) {
    if(!tree) return;

    if(!tree->left && !tree->right) {
        tree->alignment = strdup(collection->sequences[tree->id].sequence);
        return;
    }

    progressive_alignment(collection, tree->left);
    progressive_alignment(collection, tree->right);
    
    char* aligned = needleman_wunsch_align(tree->left->alignment, tree->right->alignment);
    free(tree->alignment);
    tree->alignment = aligned;
}
