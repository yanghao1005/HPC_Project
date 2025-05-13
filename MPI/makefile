CC = mpicc
CFLAGS = -O2 -std=c99 -lm
TARGET = heat_mpi
SRC = heat_mpi.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET) *.bmp
