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
	b[id] = a[id] + 1;
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

	cudaMemcpy(model.image, b, sizeof(int) * length,cudaMemcpyDeviceToHost);

	float Time_Elapse;

	printf("time= %lf s FPS= %lf\n", Time_Elapse/1000, 1000 / Time_Elapse);
}