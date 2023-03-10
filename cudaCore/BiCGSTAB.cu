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
			if(kv.second != 0)
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

__global__ void computeP(Scalar* p,Scalar* r,Scalar* v,int length,Scalar beta,Scalar w)
{
	auto id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id >= length)
		return;
	p[id] = r[id] + beta * (p[id] - w * v[id]);
}

__global__ void MatrixMultVector(Scalar* v1,Scalar* v2,IndexValue** matrix,int* preA,int length)
{
	auto id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id >= length)
		return;
	v1[id] = 0;
	for (int i = 0; i < preA[id]; i++)
	{
		v1[id] += matrix[id][i].value * v2[matrix[id][i].colid];
	}
}

__global__ void computeS(Scalar* s,Scalar* r,Scalar* v,int length)
{
	auto id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id >= length)
		return;
	s[id] = v[id];
}

__global__ void computeX(Scalar* x,Scalar* p,Scalar* s,int length,double alpha,double w)
{
	auto id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id >= length)
		return;
	x[id] = x[id] + alpha * p[id] + w * s[id];
}

__global__ void computeR(Scalar* r,Scalar* s,Scalar* t,int length,double w)
{
	auto id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id >= length)
		return;
	r[id] = s[id] - w * t[id];
}

void BiCGSTAB(const SymetrixSparseMatrix& A,Vector& x,const Vector& b,double tolerance,int limit,int& iter,double& norm)
{
	int bs = (A.get_row() / 32 * 32 == A.get_row()) ?  A.get_row() / 32 : A.get_row() / 32 + 1;
	dim3 blockSize(bs);
	dim3 threadSize(32);
	Scalar rho0,w,alpha,rho1;
	rho0 = w = alpha = 1.0;
<<<<<<< HEAD
	thrust::device_vector<Scalar> r(b.begin(),b.end()),r_hat = r,v(b.size()),p(b.size()),temp(b.size());
	thrust::host_vector<Scalar> vec3 = r;
	// for (int i = 0; i < vec3.size(); i++)--
	
	// {
	// 	std::cout << vec3[i] << " " << vec3[i] << std::endl;
	// }
	rho1 = thrust::reduce(r.begin(),r.end());
=======
	CudaSPMatrix cspm(A.get_row(),A.get_col(),A);
	thrust::device_vector<Scalar> r(b.begin(),b.end()),xx(b.size()),r_hat = r,v(b.size()),p(b.size()),s(b.size()),t(b.size()),temp(b.size());

>>>>>>> 0ebf6ca0fe0669373ec64d4fbffec72b4f299699
	iter = 0;
	norm = 1000;
	double normb = b.norm1();
	while(iter < limit && norm > tolerance * normb)
	{
		thrust::transform(thrust::device,r_hat.begin(),r_hat.end(),r.begin(),temp.begin(),thrust::multiplies<Scalar>());
		rho1 = thrust::reduce(thrust::device,temp.begin(),temp.end());
		double beta = rho1 / rho0 * alpha / w;
		rho0 = rho1;

		computeP<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&p[0]),thrust::raw_pointer_cast(&r[0]),thrust::raw_pointer_cast(&v[0]),p.size(),beta,w);
		std::vector<Scalar> temp2 (p.begin(),p.end());
		MatrixMultVector<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&v[0]),thrust::raw_pointer_cast(&p[0]),cspm.dev_matrix,cspm.preA,v.size());
		temp2 = {v.begin(),v.end()};

		thrust::transform(thrust::device,r_hat.begin(),r_hat.end(),v.begin(),temp.begin(),thrust::multiplies<Scalar>());
		temp2 = {temp.begin(),temp.end()};
		alpha = rho1 / thrust::reduce(thrust::device,temp.begin(),temp.end());

		computeS<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&s[0]),thrust::raw_pointer_cast(&r[0]),thrust::raw_pointer_cast(&v[0]),s.size());
		temp2 = {s.begin(),s.end()};
		MatrixMultVector<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&t[0]),thrust::raw_pointer_cast(&s[0]),cspm.dev_matrix,cspm.preA,t.size());
		temp2 = {t.begin(),t.end()};


		thrust::transform(thrust::device,s.begin(),s.end(),t.begin(),temp.begin(),thrust::multiplies<Scalar>());
		w = thrust::reduce(thrust::device,temp.begin(),temp.end());
		temp2 = {temp.begin(),temp.end()};
		thrust::transform(thrust::device,t.begin(),t.end(),t.begin(),temp.begin(),thrust::multiplies<Scalar>());
		temp2 = {temp.begin(),temp.end()};
		w = w / thrust::reduce(thrust::device,temp.begin(),temp.end());


		// thrust::transform(thrust::device,r_hat.begin(),r_hat.end(),t.begin(),temp.begin(),thrust::multiplies<Scalar>());
		// temp2 = {temp.begin(),temp.end()};
		// rho1 = -w * thrust::reduce(thrust::device,temp.begin(),temp.end());
		computeX<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&xx[0]),thrust::raw_pointer_cast(&p[0]),thrust::raw_pointer_cast(&s[0]),xx.size(),alpha,w);
		temp2 = {xx.begin(),xx.end()};
		computeR<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&r[0]),thrust::raw_pointer_cast(&s[0]),thrust::raw_pointer_cast(&t[0]),r.size(),w);
		temp2 = {r.begin(),r.end()};
		iter++;
		thrust::transform(thrust::device,r.begin(),r.end(),r.begin(),temp.begin(),thrust::multiplies<Scalar>());
		temp2 = {temp.begin(),temp.end()};
		norm = thrust::reduce(thrust::device,temp.begin(),temp.end());
		norm = std::sqrt(norm);
	}
	x.setvalues({xx.begin(),xx.end()});
}