#!/usr/bin/env bash
# Parallel build script for the project using OpenMPI

# Create build and bin directories if they don't exist
mkdir -p build/obj/parallel/
mkdir -p bin

# Compile each source file into an object file in the 'build' directory
for src_file in src/*.c; do
    obj_file="build/obj/parallel/$(basename "$src_file" .c).o"
    mpicc -Wall -g -Iinclude -c "$src_file" -o "$obj_file"
done

# Link the object files into the final parallel executable
mpicc -Wall -g -Iinclude -o bin/msafft_parallel build/obj/parallel/*.o

echo "Parallel build complete. Executable is in bin/msafft_parallel"
