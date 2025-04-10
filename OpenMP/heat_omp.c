#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define BMP_HEADER_SIZE 54
#define ALPHA 0.01      // Thermal diffusivity
#define L 0.2           // Length (m) of the square domain
#define DX 0.02         // Grid spacing in x-direction
#define DY 0.02         // Grid spacing in y-direction
#define DT 0.0005       // Time step
#define T 1500          // Temperature (K) for the heat source

// Optional: Function to print the grid (for debugging or visualization)
void print_grid(double *grid, int nx, int ny) {
    int i, j;
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            printf("%.2f ", grid[i * ny + j]);
        }
        printf("\n");
    }
    printf("\n");
}

// Function to initialize the grid with a heat source along both diagonals
void initialize_grid(double *grid, int nx, int ny, int temp_source) {
    int i, j;
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            if (i == j || i == nx - 1 - j)
                grid[i * ny + j] = (double) temp_source;
            else
                grid[i * ny + j] = 0.0;
        }
    }
}

// Function to solve the heat equation using finite differences and OpenMP parallelization
void solve_heat_equation(double *grid, double *new_grid, int steps, double r, int nx, int ny) {
    int step, i, j;
    for (step = 0; step < steps; step++) {
        // Update the interior points (embarrassingly parallel)
        #pragma omp parallel for collapse(2) default(none) shared(grid, new_grid, nx, ny, r)
        for (i = 1; i < nx - 1; i++) {
            for (j = 1; j < ny - 1; j++) {
                new_grid[i * ny + j] = grid[i * ny + j] +
                    r * (grid[(i + 1) * ny + j] + grid[(i - 1) * ny + j] - 2.0 * grid[i * ny + j]) +
                    r * (grid[i * ny + (j + 1)] + grid[i * ny + (j - 1)] - 2.0 * grid[i * ny + j]);
            }
        }
        // Update the boundaries (Dirichlet conditions)
        #pragma omp parallel for default(none) shared(new_grid, nx, ny)
        for (i = 0; i < nx; i++) {
            new_grid[i] = 0.0;                     // Top boundary
            new_grid[ny * (nx - 1) + i] = 0.0;       // Bottom boundary
        }
        #pragma omp parallel for default(none) shared(new_grid, nx, ny)
        for (j = 0; j < ny; j++) {
            new_grid[j * nx] = 0.0;                // Left boundary
            new_grid[(ny - 1) + j * nx] = 0.0;       // Right boundary
        }
        // Swap pointers so that new_grid becomes the grid for the next time step
        double *temp = grid;
        grid = new_grid;
        new_grid = temp;
    }
}

// Function to write the BMP file header
void write_bmp_header(FILE *file, int width, int height) {
    unsigned char header[BMP_HEADER_SIZE] = {0};
    int file_size = BMP_HEADER_SIZE + 3 * width * height;
    header[0] = 'B';
    header[1] = 'M';
    header[2] = file_size & 0xFF;
    header[3] = (file_size >> 8) & 0xFF;
    header[4] = (file_size >> 16) & 0xFF;
    header[5] = (file_size >> 24) & 0xFF;
    header[10] = BMP_HEADER_SIZE;
    header[14] = 40;  // Info header size
    header[18] = width & 0xFF;
    header[19] = (width >> 8) & 0xFF;
    header[20] = (width >> 16) & 0xFF;
    header[21] = (width >> 24) & 0xFF;
    header[22] = height & 0xFF;
    header[23] = (height >> 8) & 0xFF;
    header[24] = (height >> 16) & 0xFF;
    header[25] = (height >> 24) & 0xFF;
    header[26] = 1;   // Planes
    header[28] = 24;  // Bits per pixel

    fwrite(header, 1, BMP_HEADER_SIZE, file);
}

// Function to determine RGB color based on temperature value
void get_color(double value, unsigned char *r, unsigned char *g, unsigned char *b) {
    if (value >= 500.0) {
        *r = 255; *g = 0; *b = 0;         // Red
    } else if (value >= 100.0) {
        *r = 255; *g = 128; *b = 0;         // Orange
    } else if (value >= 50.0) {
        *r = 171; *g = 71; *b = 188;        // Lilac
    } else if (value >= 25.0) {
        *r = 255; *g = 255; *b = 0;          // Yellow
    } else if (value >= 1.0) {
        *r = 0; *g = 0; *b = 255;           // Blue
    } else if (value >= 0.1) {
        *r = 5; *g = 248; *b = 252;         // Cyan
    } else {
        *r = 255; *g = 255; *b = 255;        // White
    }
}

// Function to write the temperature grid data into a BMP file
void write_grid(FILE *file, double *grid, int nx, int ny) {
    int i, j, padding;
    // BMP files store rows in reverse order (bottom-to-top)
    for (i = nx - 1; i >= 0; i--) {
        for (j = 0; j < ny; j++) {
            unsigned char r, g, b;
            get_color(grid[i * ny + j], &r, &g, &b);
            fwrite(&b, 1, 1, file); // Write blue channel
            fwrite(&g, 1, 1, file); // Write green channel
            fwrite(&r, 1, 1, file); // Write red channel
        }
        // Row padding for 4-byte alignment
        for (padding = 0; padding < (4 - (nx * 3) % 4) % 4; padding++) {
            fputc(0, file);
        }
    }
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        printf("Command line wrong\n");
        printf("Usage: %s size steps output_file.bmp\n", argv[0]);
        return 1;
    }

    int nx, ny, steps;
    double r;

    nx = ny = atoi(argv[1]);
    steps = atoi(argv[2]);
    r = ALPHA * DT / (DX * DY);

    // Start the timer using omp_get_wtime (time starts before allocating variables)
    double start_time = omp_get_wtime();

    // Allocate memory for the grid and a temporary grid for updates
    double *grid = (double *)calloc(nx * ny, sizeof(double));
    double *new_grid = (double *)calloc(nx * ny, sizeof(double));

    // Initialize the grid
    initialize_grid(grid, nx, ny, T);

    // Solve the heat equation (the main computational kernel is parallelized)
    solve_heat_equation(grid, new_grid, steps, r, nx, ny);

    // Open the output BMP file and write the header and grid data
    FILE *file = fopen(argv[3], "wb");
    if (!file) {
        printf("Error opening the output file.\n");
        free(grid);
        free(new_grid);
        return 1;
    }
    write_bmp_header(file, nx, ny);
    write_grid(file, grid, nx, ny);
    fclose(file);

    // Free allocated memory (time measurement stops after freeing variables)
    free(grid);
    free(new_grid);

    // Stop the timer after all computations and memory free operations are complete
    double end_time = omp_get_wtime();
    printf("The Execution Time = %f seconds with a matrix size of %dx%d and %d steps\n",
           end_time - start_time, nx, nx, steps);

    return 0;
}
