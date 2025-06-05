NVCC = nvcc
NVCCFLAGS = -O3 -arch=sm_75

all: heat_cuda

heat_cuda: heat_cuda.cu
	$(NVCC) $(NVCCFLAGS) -o $@ $<

clean:
	rm -f heat_cuda *.o *.bmp

