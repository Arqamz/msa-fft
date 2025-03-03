#include "project.h"
#include "utils.h"

int main(int argc, char** argv) {
    if(argc != 2) {
        fprintf(stderr, "Usage: %s <input.fasta>\n", argv[0]);
        return 1;
    }

    Timer total_timer, read_timer, fft_timer, tree_timer, align_timer;
    
    start_timer(&total_timer, "Total Program Execution");

    // Read and encode sequences
    start_timer(&read_timer, "Read FASTA and Encode Sequences");
    SequenceCollection* collection = read_fasta(argv[1]);
    stop_timer(&read_timer);
    if(!collection) return 1;
    print_elapsed(read_timer);

    // FFT-based distance matrix
    start_timer(&fft_timer, "FFT Distance Matrix Calculation");
    double** dist_matrix = compute_fft_distance_matrix(collection);
    stop_timer(&fft_timer);
    print_elapsed(fft_timer);
    print_distance_matrix(dist_matrix, collection->count);

    // Guide tree construction
    start_timer(&tree_timer, "UPGMA Guide Tree Construction");
    UPGMA_Node* tree = build_guide_tree(dist_matrix, collection->count);
    stop_timer(&tree_timer);
    print_elapsed(tree_timer);
    print_guide_tree(tree);

    // Progressive alignment
    start_timer(&align_timer, "Progressive Alignment");
    progressive_alignment(collection, tree);
    stop_timer(&align_timer);
    print_elapsed(align_timer);

    printf("\nFinal Alignment:\n%s\n", tree->alignment);

    // Cleanup
    free_matrix(dist_matrix, collection->count);
    free_tree(tree);
    free_collection(collection);
    
    stop_timer(&total_timer);
    print_elapsed(total_timer);

    return 0;
}