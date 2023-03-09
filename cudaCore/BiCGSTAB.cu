#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include<thrust/host_vector.h>
#include<thrust/device_vector.h>
#include<thrust/copy.h>
#include<thrust/sort.h>
#include<thrust/execution_policy.h>

// #include"Math/BiCGSTABSolver.hpp"
#include "BICGSTAB.cuh"

__global__ void compute(int* a,int* b,int length)
{
	auto id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id >= length)
		return;
	// b[id] = a[id] + 1;
	// b[id] = 1;
	printf("asdas");
}

void BiCGSTAB()
{
	
	dim3 blockSize(32 ,32);
	dim3 threadSize(32, 32);
	int* a,*b;
	const int length = 4;
	cudaMalloc((void**)&a, sizeof(int) * length);
	cudaMalloc((void**)&b, sizeof(int) * length);
	
	compute << <blockSize, threadSize >> > (a,b,length);
	int bb[length];

	cudaMemcpy(bb, b, sizeof(int) * length,cudaMemcpyDeviceToHost);
	cudaFree(a);
	cudaFree(b);
}