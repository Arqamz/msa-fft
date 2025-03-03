#!/usr/bin/env bash
# Serial build script for the project

# Create build and bin directories if they don't exist
mkdir -p build/obj/serial
mkdir -p bin

# Compile each source file into an object file in the 'build' directory
for src_file in src/*.c; do
    obj_file="build/obj/serial/$(basename "$src_file" .c).o"
    gcc -Wall -g -Iinclude -c "$src_file" -o "$obj_file"
done

# Link the object files into the final executable
gcc -Wall -g -Iinclude -o bin/msafft_serial build/obj/serial/*.o

echo "Serial build complete. Executable is in bin/msafft_serial"
