all: thanos
thanos:
	/usr/local/cuda/bin/nvcc -G -O3 -I/usr/local/cuda/include -std=c++11  -gencode arch=compute_60,code=sm_60 -Xcompiler -fopenmp -lnvToolsExt -lgomp  Thanos.cu -o thanos
clean:
	rm thanos


