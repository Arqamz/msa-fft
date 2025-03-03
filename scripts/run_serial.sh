#!/usr/bin/env bash
# Serial run script for the project

# Check if the executable exists
if [ ! -f bin/msafft_serial ]; then
  echo "Error: Serial executable not found. Please build it first."
  exit 1
fi

# Run the serial executable
bin/msafft_serial

echo "Serial program run complete."
