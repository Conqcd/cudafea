#include "BICGSTAB.cuh"
#include<vector>
#include<iostream>

CudaVector::~CudaVector()
{
}

void CudaVector::AllocateData(const Vector& vec)
{
	auto v = vec.generateScalar();
	thrust::copy(v.begin(),v.end(),values.begin());
}

CudaSPVector::~CudaSPVector()
{
}

void CudaSPVector::AllocateData(const Vector& vector)
{
	std::vector<idxType> ii;
	std::vector<Scalar> vv;
	for (int i = 0; i < vector.size(); i++)
	{
		if(vector[i] != 0)
			ii.push_back(i),vv.push_back(vector[i]);
	}
	
}

CudaSPMatrix::~CudaSPMatrix()
{
	for (int i = 0; i < row; i++)
	{
		cudaFree(matrix[i]);
	}
	cudaFree(dev_matrix);
	cudaFree(preA);
	delete[] matrix;
}

void CudaSPMatrix::AllocateData(const SymetrixSparseMatrix& mat)
{
	matrix = new IndexValue*[row];
	auto prea = new int[row];
	cudaMalloc((void**)&preA, sizeof(int) * row);
	for (int i = 0; i < row; i++)
	{
		std::vector<IndexValue> vv;
		for (auto& kv:mat.getRow(i))
		{
			if(kv.second != 0 && kv.first >= row)
				vv.push_back({(int)kv.first,kv.second});
		}
		prea[i] = vv.size();
		cudaMalloc((void**)&matrix[i], sizeof(IndexValue) * vv.size());
		cudaMemcpy(matrix[i], vv.data(), sizeof(IndexValue) * vv.size(),cudaMemcpyHostToDevice);
	}
	cudaMemcpy(preA, prea, sizeof(int) * row,cudaMemcpyHostToDevice);
	delete []prea;
	cudaMalloc((void**)&dev_matrix, sizeof(IndexValue*) * row);
	cudaMemcpy(dev_matrix, matrix, sizeof(IndexValue*) * row,cudaMemcpyHostToDevice);
}

__global__ void compute(thrust::pair<int,double>* a,int length)
{
	auto id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id >= length)
		return;
	// b[id] = a[id] + 1;
	// b[id] = 1;
	a[id].first = id + 3;
	a[id].second = id + 2;
}

void BiCGSTAB(const SymetrixSparseMatrix& A,Vector& x,const Vector& b,double tolerance,int limit,int& iter,double& norm)
{
	int bs = std::sqrt(A.get_row());
	int th = (A.get_row() / bs * bs == A.get_row()) ?  A.get_row() / bs : A.get_row() / bs+ 1;
	dim3 blockSize(bs);
	dim3 threadSize(th);
	double rho0,w,alpha,rho1;
	rho0 = w = alpha = 1.0;
	thrust::device_vector<Scalar> r(b.begin(),b.end()),r_hat = r,v(b.size()),p(b.size()),temp(b.size());
	thrust::host_vector<Scalar> vec3 = r;
	// for (int i = 0; i < vec3.size(); i++)--
	
	// {
	// 	std::cout << vec3[i] << " " << vec3[i] << std::endl;
	// }
	rho1 = thrust::reduce(r.begin(),r.end());
	iter = 0;
	norm = 1000;
	double normb = b.norm1();;
	while(iter < limit && norm > tolerance * normb)
	{
		double beta = rho1 / rho0 * alpha / w;

		iter++;
	}
	// compute<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&vec[0]),vec.size());
	// vec3 = r;
	for (int i = 0; i < vec3.size(); i++)
	{
		std::cout << vec3[i] << vec3[i] << std::endl;
	}
}