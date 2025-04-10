CC = gcc
CFLAGS = -fopenmp -O2 -lm
TARGET = heat_omp
SRC = heat_omp.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET) *.bmp
