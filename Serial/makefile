# Makefile for compiling heat_serial.c with gcc

# Compiler
CC = gcc

# Compiler flags
CFLAGS = -Wall -O2

# Target executable
TARGET = heat_serial

# Source file
SRC = heat_serial.c

# Default target
all: $(TARGET)

# Build target
$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC)

# Clean up build files
clean:
	rm -f $(TARGET)