#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>

#define BMP_HEADER_SIZE 54
#define ALPHA 0.01 // Thermal diffusivity
#define L 0.2      // Length (m) of the square domain
#define DX 0.02    // Grid spacing in x-direction
#define DY 0.02    // Grid spacing in y-direction
#define DT 0.0005  // Time step
#define T 1500     // Temperature (K) for the heat source

// Optional: Function to print the grid (for debugging or visualization)
void print_grid(double *grid, int nx, int ny)
{
    int i, j;
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            printf("%.2f ", grid[i * ny + j]);
        }
        printf("\n");
    }
    printf("\n");
}

// Function to initialize a portion of the grid with a heat source along both diagonals
void initialize_grid_portion(double *grid, int start_row, int end_row, int nx, int ny, int temp_source)
{
    int i, j;
    for (i = start_row; i < end_row; i++)
    {
        for (j = 0; j < ny; j++)
        {
            if (i == j || i == nx - 1 - j)
                grid[(i - start_row) * ny + j] = (double)temp_source;
            else
                grid[(i - start_row) * ny + j] = 0.0;
        }
    }
}

// Function to solve the heat equation using finite differences with hybrid OpenMP/MPI parallelization
void solve_heat_equation_hybrid(double *grid, double *new_grid, int steps, double r,
                                int local_nx, int ny, int rank, int size, MPI_Comm comm)
{
    int step, i, j;
    // Allocate buffers for sending and receiving halo regions
    double *send_top = (double *)malloc(ny * sizeof(double));
    double *send_bottom = (double *)malloc(ny * sizeof(double));
    double *recv_top = (double *)malloc(ny * sizeof(double));
    double *recv_bottom = (double *)malloc(ny * sizeof(double));

    // MPI request handles for non-blocking communication
    MPI_Request requests[4];
    MPI_Status statuses[4];

    for (step = 0; step < steps; step++)
    {
        int req_count = 0;

        // Exchange halo regions with neighboring processes
        // Send top row to rank-1 and receive from rank-1
        if (rank > 0)
        {
            // Prepare top row to send
            for (j = 0; j < ny; j++)
            {
                send_top[j] = grid[j];
            }
            MPI_Isend(send_top, ny, MPI_DOUBLE, rank - 1, 0, comm, &requests[req_count++]);
            MPI_Irecv(recv_top, ny, MPI_DOUBLE, rank - 1, 1, comm, &requests[req_count++]);
        }

        // Send bottom row to rank+1 and receive from rank+1
        if (rank < size - 1)
        {
            // Prepare bottom row to send
            for (j = 0; j < ny; j++)
            {
                send_bottom[j] = grid[(local_nx - 1) * ny + j];
            }
            MPI_Isend(send_bottom, ny, MPI_DOUBLE, rank + 1, 1, comm, &requests[req_count++]);
            MPI_Irecv(recv_bottom, ny, MPI_DOUBLE, rank + 1, 0, comm, &requests[req_count++]);
        }

        // Wait for all communications to complete
        if (req_count > 0)
        {
            MPI_Waitall(req_count, requests, statuses);
        }
    }

// Update the interior points using OpenMP (excluding halo rows)
#pragma omp parallel for collapse(2) default(none) shared(grid, new_grid, local_nx, ny, r, rank, size)
    for (i = 1; i < local_nx - 1; i++)
    {
        for (j = 1; j < ny - 1; j++)
        {
            new_grid[i * ny + j] = grid[i * ny + j] +
                                   r * (grid[(i + 1) * ny + j] + grid[(i - 1) * ny + j] - 2.0 * grid[i * ny + j]) +
                                   r * (grid[i * ny + (j + 1)] + grid[i * ny + (j - 1)] - 2.0 * grid[i * ny + j]);
        }
    }

    // Update the boundary points that depend on halo regions
    if (rank > 0)
    {
// Update top row (i=0) using data received from rank-1
#pragma omp parallel for default(none) shared(grid, new_grid, ny, r, recv_top)
        for (j = 1; j < ny - 1; j++)
        {
            new_grid[j] = grid[j] +
                          r * (grid[ny + j] + recv_top[j] - 2.0 * grid[j]) +
                          r * (grid[j + 1] + grid[j - 1] - 2.0 * grid[j]);
        }
    }
    else
    {
// Process 0 has the top global boundary (Dirichlet condition = 0)
#pragma omp parallel for default(none) shared(new_grid, ny)
        for (j = 0; j < ny; j++)
        {
            new_grid[j] = 0.0;
        }
    }

    if (rank < size - 1)
    {
// Update bottom row (i=local_nx-1) using data received from rank+1
#pragma omp parallel for default(none) shared(grid, new_grid, local_nx, ny, r, recv_bottom)
        for (j = 1; j < ny - 1; j++)
        {
            new_grid[(local_nx - 1) * ny + j] = grid[(local_nx - 1) * ny + j] +
                                                r * (recv_bottom[j] + grid[(local_nx - 2) * ny + j] - 2.0 * grid[(local_nx - 1) * ny + j]) +
                                                r * (grid[(local_nx - 1) * ny + (j + 1)] + grid[(local_nx - 1) * ny + (j - 1)] - 2.0 * grid[(local_nx - 1) * ny + j]);
        }
    }
    else
    {
// Last process has the bottom global boundary (Dirichlet condition = 0)
#pragma omp parallel for default(none) shared(new_grid, local_nx, ny)
        for (j = 0; j < ny; j++)
        {
            new_grid[(local_nx - 1) * ny + j] = 0.0;
        }
    }

// Update left and right boundaries (Dirichlet conditions = 0)
#pragma omp parallel for default(none) shared(new_grid, local_nx, ny)
    for (i = 0; i < local_nx; i++)
    {
        new_grid[i * ny] = 0.0;            // Left boundary
        new_grid[i * ny + (ny - 1)] = 0.0; // Right boundary
    }

    // Swap pointers so that new_grid becomes the grid for the next time step
    double *temp = grid;
    grid = new_grid;
    new_grid = temp;

    // Free allocated buffers
    free(send_top);
    free(send_bottom);
    free(recv_top);
    free(recv_bottom);
}

// Function to write the BMP file header
void write_bmp_header(FILE *file, int width, int height)
{
    unsigned char header[BMP_HEADER_SIZE] = {0};
    int file_size = BMP_HEADER_SIZE + 3 * width * height;
    header[0] = 'B';
    header[1] = 'M';
    header[2] = file_size & 0xFF;
    header[3] = (file_size >> 8) & 0xFF;
    header[4] = (file_size >> 16) & 0xFF;
    header[5] = (file_size >> 24) & 0xFF;
    header[10] = BMP_HEADER_SIZE;
    header[14] = 40; // Info header size
    header[18] = width & 0xFF;
    header[19] = (width >> 8) & 0xFF;
    header[20] = (width >> 16) & 0xFF;
    header[21] = (width >> 24) & 0xFF;
    header[22] = height & 0xFF;
    header[23] = (height >> 8) & 0xFF;
    header[24] = (height >> 16) & 0xFF;
    header[25] = (height >> 24) & 0xFF;
    header[26] = 1;  // Planes
    header[28] = 24; // Bits per pixel

    fwrite(header, 1, BMP_HEADER_SIZE, file);
}

// Function to determine RGB color based on temperature value
void get_color(double value, unsigned char *r, unsigned char *g, unsigned char *b)
{
    if (value >= 500.0)
    {
        *r = 255;
        *g = 0;
        *b = 0; // Red
    }
    else if (value >= 100.0)
    {
        *r = 255;
        *g = 128;
        *b = 0; // Orange
    }
    else if (value >= 50.0)
    {
        *r = 171;
        *g = 71;
        *b = 188; // Lilac
    }
    else if (value >= 25.0)
    {
        *r = 255;
        *g = 255;
        *b = 0; // Yellow
    }
    else if (value >= 1.0)
    {
        *r = 0;
        *g = 0;
        *b = 255; // Blue
    }
    else if (value >= 0.1)
    {
        *r = 5;
        *g = 248;
        *b = 252; // Cyan
    }
    else
    {
        *r = 255;
        *g = 255;
        *b = 255; // White
    }
}

// Function to write the temperature grid data into a BMP file
void write_grid(FILE *file, double *grid, int nx, int ny)
{
    int i, j, padding;
    // BMP files store rows in reverse order (bottom-to-top)
    for (i = nx - 1; i >= 0; i--)
    {
        for (j = 0; j < ny; j++)
        {
            unsigned char r, g, b;
            get_color(grid[i * ny + j], &r, &g, &b);
            fwrite(&b, 1, 1, file); // Write blue channel
            fwrite(&g, 1, 1, file); // Write green channel
            fwrite(&r, 1, 1, file); // Write red channel
        }
        // Row padding for 4-byte alignment
        for (padding = 0; padding < (4 - (ny * 3) % 4) % 4; padding++)
        {
            fputc(0, file);
        }
    }
}

int main(int argc, char *argv[])
{
    int rank, size, provided;

    // Initialize MPI with thread support for OpenMP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Check command-line arguments
    if (argc != 4)
    {
        if (rank == 0)
        {
            printf("Command line wrong\n");
            printf("Usage: %s size steps output_file.bmp\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    int nx, ny, steps;
    double r;

    nx = ny = atoi(argv[1]);
    steps = atoi(argv[2]);
    r = ALPHA * DT / (DX * DY);

    // Calculate local grid dimensions (row decomposition)
    int local_nx = nx / size;
    // Handle remainder rows
    if (rank < nx % size)
    {
        local_nx++;
    }

    // Calculate starting row for this process
    int start_row = rank * (nx / size);
    if (rank < nx % size)
    {
        start_row += rank;
    }
    else
    {
        start_row += nx % size;
    }
    int end_row = start_row + local_nx;

    // Start the timer
    double start_time;
    if (rank == 0)
    {
        start_time = MPI_Wtime();
    }

    // Allocate memory for the local grid and a temporary grid for updates
    double *local_grid = (double *)calloc(local_nx * ny, sizeof(double));
    double *local_new_grid = (double *)calloc(local_nx * ny, sizeof(double));

    // Initialize the local portion of the grid
    initialize_grid_portion(local_grid, start_row, end_row, nx, ny, T);

    // Solve the heat equation with hybrid OpenMP/MPI parallelization
    solve_heat_equation_hybrid(local_grid, local_new_grid, steps, r, local_nx, ny, rank, size, MPI_COMM_WORLD);

    // Gather results from all processes to the root process (rank 0)
    double *global_grid = NULL;
    if (rank == 0)
    {
        global_grid = (double *)calloc(nx * ny, sizeof(double));
    }

    // First determine receive counts and displacements for gatherv
    int *recvcounts = NULL;
    int *displs = NULL;
    if (rank == 0)
    {
        recvcounts = (int *)malloc(size * sizeof(int));
        displs = (int *)malloc(size * sizeof(int));

        int disp = 0;
        for (int i = 0; i < size; i++)
        {
            int rows_for_proc = nx / size;
            if (i < nx % size)
            {
                rows_for_proc++;
            }
            recvcounts[i] = rows_for_proc * ny;
            displs[i] = disp;
            disp += recvcounts[i];
        }
    }

    // Gather local grids to the global grid on rank 0
    MPI_Gatherv(local_grid, local_nx * ny, MPI_DOUBLE,
                global_grid, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Write output to BMP file (only on rank 0)
    if (rank == 0)
    {
        FILE *file = fopen(argv[3], "wb");
        if (!file)
        {
            printf("Error opening the output file.\n");
            free(global_grid);
            free(recvcounts);
            free(displs);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        write_bmp_header(file, ny, nx);
        write_grid(file, global_grid, nx, ny);
        fclose(file);

        // Stop the timer and print execution time
        double end_time = MPI_Wtime();
        printf("The Execution Time = %f seconds with a matrix size of %dx%d and %d steps\n",
               end_time - start_time, nx, nx, steps);

        // Free additional resources
        free(global_grid);
        free(recvcounts);
        free(displs);
    }

    // Free allocated memory
    free(local_grid);
    free(local_new_grid);

    // Finalize MPI
    MPI_Finalize();
    return 0;
}