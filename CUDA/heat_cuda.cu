#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda_runtime.h>

#define BMP_HEADER_SIZE 54
#define ALPHA 0.01 // Thermal diffusivity
#define L 0.2      // Length (m) of the square domain
#define DX 0.02    // Grid spacing in x-direction
#define DY 0.02    // Grid spacing in y-direction
#define DT 0.0005  // Time step
#define T 1500     // Temperature (K) for the heat source

#define BLOCK_SIZE 16 // Block size for CUDA kernels

// Macro for CUDA error checking
#define CUDA_CHECK(call)                                                                       \
    do                                                                                         \
    {                                                                                          \
        cudaError_t err = call;                                                                \
        if (err != cudaSuccess)                                                                \
        {                                                                                      \
            printf("CUDA error at %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
            exit(1);                                                                           \
        }                                                                                      \
    } while (0)

// CUDA kernel to initialize the grid with heat sources along both diagonals
__global__ void initialize_grid_kernel(double *grid, int nx, int ny, double temp_source)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < nx && j < ny)
    {
        int idx = i * ny + j;
        if (i == j || i == nx - 1 - j)
        {
            grid[idx] = temp_source;
        }
        else
        {
            grid[idx] = 0.0;
        }
    }
}

// CUDA kernel to update interior points of the grid
__global__ void update_interior_kernel(double *grid, double *new_grid, int nx, int ny, double r)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    // Only update interior points (not boundaries)
    if (i > 0 && i < nx - 1 && j > 0 && j < ny - 1)
    {
        int idx = i * ny + j;
        new_grid[idx] = grid[idx] +
                        r * (grid[(i + 1) * ny + j] + grid[(i - 1) * ny + j] - 2.0 * grid[idx]) +
                        r * (grid[i * ny + (j + 1)] + grid[i * ny + (j - 1)] - 2.0 * grid[idx]);
    }
}

// CUDA kernel to apply boundary conditions (Dirichlet: u=0 on boundaries)
__global__ void apply_boundaries_kernel(double *new_grid, int nx, int ny)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < nx)
    {
        // Top and bottom boundaries
        new_grid[idx] = 0.0;                 // Top boundary
        new_grid[(nx - 1) * ny + idx] = 0.0; // Bottom boundary
    }

    if (idx < ny)
    {
        // Left and right boundaries
        new_grid[idx * nx] = 0.0;            // Left boundary
        new_grid[idx * nx + (ny - 1)] = 0.0; // Right boundary
    }
}

// Function to solve the heat equation using CUDA
void solve_heat_equation_cuda(double *h_grid, int steps, double r, int nx, int ny, int block_x, int block_y)
{
    size_t grid_size = nx * ny * sizeof(double);

    // Allocate device memory
    double *d_grid, *d_new_grid;
    CUDA_CHECK(cudaMalloc(&d_grid, grid_size));
    CUDA_CHECK(cudaMalloc(&d_new_grid, grid_size));

    // Copy initial grid to device
    CUDA_CHECK(cudaMemcpy(d_grid, h_grid, grid_size, cudaMemcpyHostToDevice));

    // Define grid and block dimensions with custom block size
    dim3 blockSize(block_x, block_y);
    dim3 gridSize((nx + blockSize.x - 1) / blockSize.x,
                  (ny + blockSize.y - 1) / blockSize.y);

    // Linear grid for boundary kernel
    dim3 linearBlockSize(256);
    dim3 linearGridSize((max(nx, ny) + linearBlockSize.x - 1) / linearBlockSize.x);

    printf("Using block size: %dx%d (%d threads per block)\n", block_x, block_y, block_x * block_y);
    printf("Grid configuration: %dx%d blocks\n", gridSize.x, gridSize.y);
    printf("Total threads: %d\n", gridSize.x * gridSize.y * block_x * block_y);

    // Time stepping loop
    for (int step = 0; step < steps; step++)
    {
        // Update interior points
        update_interior_kernel<<<gridSize, blockSize>>>(d_grid, d_new_grid, nx, ny, r);
        CUDA_CHECK(cudaGetLastError());

        // Apply boundary conditions
        apply_boundaries_kernel<<<linearGridSize, linearBlockSize>>>(d_new_grid, nx, ny);
        CUDA_CHECK(cudaGetLastError());

        // Synchronize every 50 steps instead of every step for better performance
        if (step % 50 == 0 || step == steps - 1)
        {
            CUDA_CHECK(cudaDeviceSynchronize());
        }

        // Swap pointers
        double *temp = d_grid;
        d_grid = d_new_grid;
        d_new_grid = temp;
    }

    // Final synchronization
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy result back to host
    CUDA_CHECK(cudaMemcpy(h_grid, d_grid, grid_size, cudaMemcpyDeviceToHost));

    // Free device memory
    CUDA_CHECK(cudaFree(d_grid));
    CUDA_CHECK(cudaFree(d_new_grid));
}

// Function to initialize the grid on CPU (for verification)
void initialize_grid_cpu(double *grid, int nx, int ny, double temp_source)
{
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            if (i == j || i == nx - 1 - j)
            {
                grid[i * ny + j] = temp_source;
            }
            else
            {
                grid[i * ny + j] = 0.0;
            }
        }
    }
}

// Alternative GPU-only initialization
void initialize_grid_gpu(double *h_grid, int nx, int ny, double temp_source)
{
    size_t grid_size = nx * ny * sizeof(double);

    // Allocate device memory
    double *d_grid;
    CUDA_CHECK(cudaMalloc(&d_grid, grid_size));

    // Define grid and block dimensions
    dim3 blockSize(BLOCK_SIZE, BLOCK_SIZE);
    dim3 gridSize((nx + blockSize.x - 1) / blockSize.x,
                  (ny + blockSize.y - 1) / blockSize.y);

    // Initialize grid on GPU
    initialize_grid_kernel<<<gridSize, blockSize>>>(d_grid, nx, ny, temp_source);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy result back to host
    CUDA_CHECK(cudaMemcpy(h_grid, d_grid, grid_size, cudaMemcpyDeviceToHost));

    // Free device memory
    CUDA_CHECK(cudaFree(d_grid));
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

void print_grid(double *grid, int nx, int ny)
{
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            printf("%.2f ", grid[i * ny + j]);
        }
        printf("\n");
    }
    printf("\n");
}

int main(int argc, char *argv[])
{
    // Check command-line arguments
    if (argc != 6)
    {
        printf("Command line wrong\n");
        printf("Usage: %s size steps block_x block_y output_file.bmp\n", argv[0]);
        printf("Example: %s 512 1000 16 16 output.bmp\n", argv[0]);
        return 1;
    }

    int nx, ny, steps, block_x, block_y;
    double r;

    nx = ny = atoi(argv[1]);
    steps = atoi(argv[2]);
    block_x = atoi(argv[3]);
    block_y = atoi(argv[4]);
    r = ALPHA * DT / (DX * DY);

    // Validate block size
    if (block_x * block_y > 1024)
    {
        printf("Error: Block size %dx%d = %d exceeds maximum 1024 threads per block\n",
               block_x, block_y, block_x * block_y);
        return 1;
    }

    printf("========== CUDA HEAT EQUATION SOLVER ==========\n");
    printf("Grid size: %dx%d\n", nx, ny);
    printf("Time steps: %d\n", steps);
    printf("Block size: %dx%d (%d threads)\n", block_x, block_y, block_x * block_y);
    printf("Output file: %s\n", argv[5]);
    printf("===============================================\n");

    // Check CUDA device availability
    int deviceCount;
    CUDA_CHECK(cudaGetDeviceCount(&deviceCount));
    if (deviceCount == 0)
    {
        printf("No CUDA devices found!\n");
        return 1;
    }

    // Print device information
    cudaDeviceProp deviceProp;
    CUDA_CHECK(cudaGetDeviceProperties(&deviceProp, 0));
    printf("Using GPU: %s\n", deviceProp.name);
    printf("Max threads per block: %d\n", deviceProp.maxThreadsPerBlock);
    printf("Number of SMs: %d\n", deviceProp.multiProcessorCount);

    // START TIMING HERE - Same position as Serial/OpenMP versions
    cudaEvent_t start, stop;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));
    CUDA_CHECK(cudaEventRecord(start));

    // Allocate memory for the grid
    double *grid = (double *)malloc(nx * ny * sizeof(double));
    if (!grid)
    {
        printf("Failed to allocate host memory!\n");
        return 1;
    }

    // Initialize the grid
    initialize_grid_cpu(grid, nx, ny, T);

    // Solve the heat equation using CUDA with custom block size
    solve_heat_equation_cuda(grid, steps, r, nx, ny, block_x, block_y);

    // Write output to BMP file
    FILE *file = fopen(argv[5], "wb");
    if (!file)
    {
        printf("Error opening the output file.\n");
        free(grid);
        return 1;
    }

    write_bmp_header(file, ny, nx);
    write_grid(file, grid, nx, ny);
    fclose(file);

    // Free allocated memory
    free(grid);

    // STOP TIMING HERE - Same position as Serial/OpenMP versions
    CUDA_CHECK(cudaEventRecord(stop));
    CUDA_CHECK(cudaEventSynchronize(stop));

    float milliseconds = 0;
    CUDA_CHECK(cudaEventElapsedTime(&milliseconds, start, stop));

    // Calculate performance metrics
    long long total_flops = (long long)steps * (nx - 2) * (ny - 2) * 7; // 7 FLOPs per interior point
    double gflops = (total_flops / 1e9) / (milliseconds / 1000.0);

    long long total_bytes = (long long)steps * nx * ny * sizeof(double) * 6; // 6 memory accesses per point
    double bandwidth = (total_bytes / 1e9) / (milliseconds / 1000.0);

    // Print detailed results
    printf("\n========== PERFORMANCE RESULTS ==========\n");
    printf("Execution Time: %.4f seconds\n", milliseconds / 1000.0);
    printf("Performance: %.2f GFLOP/s\n", gflops);
    printf("Memory Bandwidth: %.2f GB/s\n", bandwidth);
    printf("Grid blocks: %dx%d = %d total blocks\n",
           (nx + block_x - 1) / block_x, (ny + block_y - 1) / block_y,
           ((nx + block_x - 1) / block_x) * ((ny + block_y - 1) / block_y));
    printf("=========================================\n");

    // Print execution time in same format as other versions
    printf("The Execution Time = %f seconds with a matrix size of %dx%d and %d steps\n",
           milliseconds / 1000.0, nx, nx, steps);

    // Cleanup CUDA events
    CUDA_CHECK(cudaEventDestroy(start));
    CUDA_CHECK(cudaEventDestroy(stop));

    return 0;
}