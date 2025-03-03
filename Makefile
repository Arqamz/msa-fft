PROJECT_NAME_SERIAL = msafft_serial
PROJECT_NAME_PARALLEL = msafft_parallel

# Set the directories
BUILD_DIR = build
BIN_DIR = bin
SRC_DIR = src
INCLUDE_DIR = include
OBJ_DIR = $(BUILD_DIR)/obj

# Set the compiler and flags
CC = gcc
MPI_CC = mpicc  # OpenMPI compiler for parallel compilation
CFLAGS = -Wall -g -I$(INCLUDE_DIR)

# Find all the source files
SOURCES = $(wildcard $(SRC_DIR)/*.c)

# Separate object files for serial and parallel builds
SERIAL_OBJECTS = $(SOURCES:$(SRC_DIR)/%.c=$(OBJ_DIR)/serial/%.o)
PARALLEL_OBJECTS = $(SOURCES:$(SRC_DIR)/%.c=$(OBJ_DIR)/parallel/%.o)

# Default target (build serial)
all: build_serial

# Create the bin directory if it doesn't exist
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Create the build directory if it doesn't exist
$(OBJ_DIR)/serial:
	mkdir -p $(OBJ_DIR)/serial

$(OBJ_DIR)/parallel:
	mkdir -p $(OBJ_DIR)/parallel

# Build serial: compile and link the serial executable
build_serial: $(BIN_DIR)/$(PROJECT_NAME_SERIAL)

# Rule to link the object files to create the serial executable
$(BIN_DIR)/$(PROJECT_NAME_SERIAL): $(SERIAL_OBJECTS) | $(BIN_DIR)
	$(CC) $(SERIAL_OBJECTS) -o $@

# Rule to compile source files into object files (serial)
$(OBJ_DIR)/serial/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)/serial
	$(CC) $(CFLAGS) -c $< -o $@

# Build parallel: compile and link the parallel executable with OpenMPI
build_parallel: $(BIN_DIR)/$(PROJECT_NAME_PARALLEL)

# Rule to link the object files to create the parallel executable
$(BIN_DIR)/$(PROJECT_NAME_PARALLEL): $(PARALLEL_OBJECTS) | $(BIN_DIR)
	$(MPI_CC) $(PARALLEL_OBJECTS) -o $@

# Rule to compile source files into object files (parallel)
$(OBJ_DIR)/parallel/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)/parallel
	$(MPI_CC) $(CFLAGS) -c $< -o $@

# Clean up all generated files
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

# Run the program serially
run_serial: $(BIN_DIR)/$(PROJECT_NAME_SERIAL)
	$(BIN_DIR)/$(PROJECT_NAME_SERIAL)

# Run the program in parallel using OpenMPI
# make run_parallel NP=4 MPI_FLAGS="--bind-to none --use-hwthread-cpus --hostfile myhostfile"
run_parallel: $(BIN_DIR)/$(PROJECT_NAME_PARALLEL)
	mpirun -np $(NP) $(MPI_FLAGS) $(BIN_DIR)/$(PROJECT_NAME_PARALLEL)
