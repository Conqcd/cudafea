#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include<thrust/host_vector.h>
#include<thrust/device_vector.h>
#include<thrust/copy.h>
#include<thrust/sort.h>
#include<thrust/execution_policy.h>

#include<vector>
#include "BICGSTAB.cuh"

CudaVector::~CudaVector()
{
	cudaFree(values);
}

void CudaVector::AllocateData(const Vector& vec)
{
	cudaMalloc((void**)&values, sizeof(double) * row);
	cudaMemcpy(values, vec.generateScalar().data(), sizeof(double) * vec.size(),cudaMemcpyHostToDevice);
}

CudaSPVector::~CudaSPVector()
{
	cudaFree(vec);
}

void CudaSPVector::AllocateData(const Vector& vector)
{
	std::vector<IndexValue> vv;
	for (int i = 0; i < vector.size(); i++)
	{
		if(vector[i] != 0)
			vv.push_back({i,vector[i]});
	}
	
	cudaMalloc((void**)&vec, sizeof(IndexValue) * vv.size());
	cudaMemcpy(vec, vv.data(), sizeof(IndexValue) * vv.size(),cudaMemcpyHostToDevice);
}

CudaSPMatrix::~CudaSPMatrix()
{
	for (int i = 0; i < row; i++)
	{
		cudaFree(matrix[i]);
	}
	cudaFree(dev_matrix);
}

void CudaSPMatrix::AllocateData(const SymetrixSparseMatrix& mat)
{
	matrix = new IndexValue*[row];
	for (int i = 0; i < row; i++)
	{
		std::vector<IndexValue> vv;
		for (auto& kv:mat.getRow(i))
		{
			if(kv.second != 0 && kv.first >= row)
				vv.push_back({(int)kv.first,kv.second});
		}
		cudaMalloc((void**)&matrix[i], sizeof(IndexValue) * vv.size());
		cudaMemcpy(matrix[i], vv.data(), sizeof(IndexValue) * vv.size(),cudaMemcpyHostToDevice);
	}
	
	cudaMalloc((void**)&dev_matrix, sizeof(IndexValue*) * row);
	cudaMemcpy(dev_matrix, matrix, sizeof(IndexValue*) * row,cudaMemcpyHostToDevice);
}

__global__ void compute(int* a,int* b,int length)
{
	auto id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id >= length)
		return;
	// b[id] = a[id] + 1;
	// b[id] = 1;
	printf("asdas");
}

void BiCGSTAB(const SymetrixSparseMatrix& A,Vector& x,const Vector& b)
{
	
	dim3 blockSize(32 ,32);
	dim3 threadSize(32, 32);
	int* a,*b2;
	const int length = 4;
	cudaMalloc((void**)&a, sizeof(int) * length);
	cudaMalloc((void**)&b2, sizeof(int) * length);
	
	compute << <blockSize, threadSize >> > (a,b2,length);
	int bb[length];

	cudaMemcpy(bb, b2, sizeof(int) * length,cudaMemcpyDeviceToHost);
	cudaFree(a);
	cudaFree(b2);
}