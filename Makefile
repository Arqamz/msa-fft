PROJECT_NAME_SERIAL = msafft_serial
PROJECT_NAME_PARALLEL = msafft_parallel

# Directory structure
BUILD_DIR = build
BIN_DIR = bin
SRC_DIR = src
INCLUDE_DIR = include
OBJ_DIR = $(BUILD_DIR)/obj

# Compiler settings
CC = gcc
MPI_CC = mpicc
CFLAGS_SERIAL = -Wall -g -I$(INCLUDE_DIR)/shared -I$(INCLUDE_DIR)/serial -lm
CFLAGS_PARALLEL = -Wall -g -I$(INCLUDE_DIR)/shared -I$(INCLUDE_DIR)/parallel -lm

# Source files
SERIAL_SOURCES = $(wildcard $(SRC_DIR)/serial/*.c) $(wildcard $(SRC_DIR)/shared/*.c)
PARALLEL_SOURCES = $(wildcard $(SRC_DIR)/parallel/*.c) $(wildcard $(SRC_DIR)/shared/*.c)

# Object files (correct pattern substitution)
SERIAL_OBJECTS = $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/serial/%.o,$(basename $(SERIAL_SOURCES)))
PARALLEL_OBJECTS = $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/parallel/%.o,$(basename $(PARALLEL_SOURCES)))

# Default target
all: build_serial

# Directory creation
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(OBJ_DIR)/serial:
	mkdir -p $(OBJ_DIR)/serial

$(OBJ_DIR)/parallel:
	mkdir -p $(OBJ_DIR)/parallel

# Compilation rules for serial
$(OBJ_DIR)/serial/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)/serial
	@mkdir -p $(@D)
	$(CC) $(CFLAGS_SERIAL) -c $< -o $@

# Compilation rules for parallel
$(OBJ_DIR)/parallel/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)/parallel
	@mkdir -p $(@D)
	$(MPI_CC) $(CFLAGS_PARALLEL) -c $< -o $@

# Linking targets
build_serial: $(BIN_DIR)/$(PROJECT_NAME_SERIAL)

$(BIN_DIR)/$(PROJECT_NAME_SERIAL): $(SERIAL_OBJECTS) | $(BIN_DIR)
	$(CC) $^ -o $@ -lm

build_parallel: $(BIN_DIR)/$(PROJECT_NAME_PARALLEL)

$(BIN_DIR)/$(PROJECT_NAME_PARALLEL): $(PARALLEL_OBJECTS) | $(BIN_DIR)
	$(MPI_CC) $^ -o $@ -lm

# Clean
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

# Run targets
run_serial: $(BIN_DIR)/$(PROJECT_NAME_SERIAL)
	./$<

run_parallel: $(BIN_DIR)/$(PROJECT_NAME_PARALLEL)
	mpirun -np $(NP) ./$<

# Phony targets
.PHONY: all clean build_serial build_parallel run_serial run_parallel serial parallel

# Shortcut targets
serial: build_serial
parallel: build_parallel