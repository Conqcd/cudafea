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

__global__ void computeP_CG(Scalar* p,Scalar* r,int length,Scalar beta)
{
	auto id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id >= length)
		return;
	p[id] = r[id] + beta * p[id];
}

__global__ void computeP_PCG(Scalar* p,Scalar* z,int length,Scalar beta)
{
	auto id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id >= length)
		return;
	p[id] = z[id] + beta * p[id];
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

__global__ void computeS(Scalar* s,Scalar* r,Scalar* v,int length,double alpha)
{
	auto id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id >= length)
		return;
	s[id] = r[id] - alpha * v[id];
}

__global__ void computeX(Scalar* x,Scalar* p,Scalar* s,int length,double alpha,double w)
{
	auto id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id >= length)
		return;
	x[id] = x[id] + alpha * p[id] + w * s[id];
}

__global__ void computeX_CG(Scalar* x,Scalar* p,int length,double alpha)
{
	auto id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id >= length)
		return;
	x[id] = x[id] + alpha * p[id] ;
}

__global__ void computeX_PCG(Scalar* x,Scalar* p,int length,double alpha)
{
	auto id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id >= length)
		return;
	x[id] = x[id] + alpha * p[id] ;
}


__global__ void computeR(Scalar* r,Scalar* s,Scalar* t,int length,double w)
{
	auto id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id >= length)
		return;
	r[id] = s[id] - w * t[id];
}

__global__ void computeR_CG(Scalar* r,Scalar* Ap,int length,double alpha)
{
	auto id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id >= length)
		return;
	r[id] = r[id] - alpha * Ap[id];
}

__global__ void computeR_PCG(Scalar* r,Scalar* w,int length,double alpha)
{
	auto id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id >= length)
		return;
	r[id] = r[id] - alpha * w[id];
}

void CG(const SymetrixSparseMatrix& A,Vector& x,const Vector& b,double tolerance,int limit,int& iter,double& norm)
{
	int bs = (A.get_row() / 32 * 32 == A.get_row()) ?  A.get_row() / 32 : A.get_row() / 32 + 1;
	dim3 blockSize(bs);
	dim3 threadSize(32);
	Scalar alpha = 0.0,rr0,rr1,beta = 0.0;

	CudaSPMatrix cspm(A.get_row(),A.get_col(),A);

	thrust::device_vector<Scalar> r(b.begin(),b.end()),xx(b.size()),p = r,Ap(b.size()),temp(b.size());

	iter = 0;
	norm = 1000;
	double normb = b.norm1();
	thrust::transform(thrust::device,r.begin(),r.end(),r.begin(),temp.begin(),thrust::multiplies<Scalar>());
	rr1 = thrust::reduce(thrust::device,temp.begin(),temp.end());


	while(iter < limit && norm > tolerance * normb)
	{
		MatrixMultVector<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&Ap[0]),thrust::raw_pointer_cast(&p[0]),cspm.dev_matrix,cspm.preA,Ap.size());
		thrust::transform(thrust::device,p.begin(),p.end(),Ap.begin(),temp.begin(),thrust::multiplies<Scalar>());
		alpha = rr1 / thrust::reduce(thrust::device,temp.begin(),temp.end());


		computeX_CG<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&xx[0]),thrust::raw_pointer_cast(&p[0]),xx.size(),alpha);

		computeR_CG<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&r[0]),thrust::raw_pointer_cast(&Ap[0]),r.size(),alpha);

		rr0 = rr1;
		thrust::transform(thrust::device,r.begin(),r.end(),r.begin(),temp.begin(),thrust::multiplies<Scalar>());
		rr1 = thrust::reduce(thrust::device,temp.begin(),temp.end());
		beta = rr1 / rr0;
		computeP_CG<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&p[0]),thrust::raw_pointer_cast(&r[0]),p.size(),beta);
		iter++;
		thrust::transform(thrust::device,r.begin(),r.end(),r.begin(),temp.begin(),thrust::multiplies<Scalar>());
		norm = thrust::reduce(thrust::device,temp.begin(),temp.end());
		norm = std::sqrt(norm);
	}
	x.setvalues({xx.begin(),xx.end()});
}

void PCG(const SymetrixSparseMatrix& A,Vector& x,const Vector& b,double tolerance,int limit,int& iter,double& norm)
{
	int bs = (A.get_row() / 32 * 32 == A.get_row()) ?  A.get_row() / 32 : A.get_row() / 32 + 1;
	dim3 blockSize(bs);
	dim3 threadSize(32);
	Scalar alpha = 0.0,rr0,rr1,beta = 0.0;

	// auto precon = A.ichol().inverse_lowertri();
	auto precon = A.ichol();
	auto preconT = precon.transpose();
	CudaSPMatrix cspm(A.get_row(),A.get_col(),A);
	CudaSPMatrix prec(precon.get_row(),precon.get_col(),precon);
	CudaSPMatrix precT(preconT.get_row(),preconT.get_col(),preconT);

	thrust::device_vector<Scalar> r(b.begin(),b.end()),xx(b.size()),p(b.size()),Ap(b.size()),temp(b.size()),
	y(b.size()),z(b.size()),w(b.size()),lastr(b.size());

	MatrixMultVector<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&y[0]),thrust::raw_pointer_cast(&r[0]),prec.dev_matrix,prec.preA,y.size());
	MatrixMultVector<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&z[0]),thrust::raw_pointer_cast(&y[0]),precT.dev_matrix,precT.preA,z.size());
	p = z;
	std::vector<Scalar> tempp{z.begin(),z.end()};
	std::vector<Scalar> tempp2{y.begin(),y.end()};
	iter = 0;
	thrust::transform(thrust::device,r.begin(),r.end(),r.begin(),temp.begin(),thrust::multiplies<Scalar>());
	norm = 1000;
	double normb = b.norm1();

	thrust::transform(thrust::device,r.begin(),r.end(),z.begin(),temp.begin(),thrust::multiplies<Scalar>());
	rr0 = thrust::reduce(thrust::device,temp.begin(),temp.end());

	while(iter < limit && norm > tolerance * normb)
	{
		MatrixMultVector<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&w[0]),thrust::raw_pointer_cast(&p[0]),cspm.dev_matrix,cspm.preA,Ap.size());

		std::vector<Scalar> tempp{w.begin(),w.end()};
		thrust::transform(thrust::device,p.begin(),p.end(),w.begin(),temp.begin(),thrust::multiplies<Scalar>());
		alpha = rr0 / thrust::reduce(thrust::device,temp.begin(),temp.end());

		computeX_PCG<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&xx[0]),thrust::raw_pointer_cast(&p[0]),xx.size(),alpha);
		tempp = {xx.begin(),xx.end()};

		lastr = r;
		computeR_PCG<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&r[0]),thrust::raw_pointer_cast(&w[0]),r.size(),alpha);
		tempp = {r.begin(),r.end()};

		MatrixMultVector<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&y[0]),thrust::raw_pointer_cast(&r[0]),prec.dev_matrix,prec.preA,y.size());
		MatrixMultVector<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&z[0]),thrust::raw_pointer_cast(&y[0]),precT.dev_matrix,precT.preA,z.size());
		tempp = {z.begin(),z.end()};


		thrust::transform(thrust::device,r.begin(),r.end(),lastr.begin(),temp.begin(),thrust::minus<Scalar>());
		thrust::transform(thrust::device,temp.begin(),temp.end(),z.begin(),temp.begin(),thrust::multiplies<Scalar>());
		rr1 = thrust::reduce(thrust::device,temp.begin(),temp.end());

		beta = rr1 / rr0;
		thrust::transform(thrust::device,r.begin(),r.end(),z.begin(),temp.begin(),thrust::multiplies<Scalar>());
		rr0 = thrust::reduce(thrust::device,temp.begin(),temp.end());
		// rr0 = rr1;

		computeP_PCG<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&p[0]),thrust::raw_pointer_cast(&z[0]),p.size(),beta);
		tempp = {p.begin(),p.end()};

		iter++;
		thrust::transform(thrust::device,r.begin(),r.end(),r.begin(),temp.begin(),thrust::multiplies<Scalar>());
		norm = thrust::reduce(thrust::device,temp.begin(),temp.end());
		norm = std::sqrt(norm);
	}
	x.setvalues({xx.begin(),xx.end()});
}

void PCG2(const SymetrixSparseMatrix& A,Vector& x,const Vector& b,double tolerance,int limit,int& iter,double& norm)
{
	int bs = (A.get_row() / 32 * 32 == A.get_row()) ?  A.get_row() / 32 : A.get_row() / 32 + 1;
	dim3 blockSize(bs);
	dim3 threadSize(32);
	Scalar alpha = 0.0,rr0,rr1,beta = 0.0;

	auto precon = A.ichol().inverse_lowertri();
	CudaSPMatrix cspm(A.get_row(),A.get_col(),A);
	CudaSPMatrix prec(precon.get_row(),precon.get_col(),precon);

	thrust::device_vector<Scalar> r(b.begin(),b.end()),xx(b.size()),p(b.size()),Ap(b.size()),temp(b.size()),z(b.size()),w(b.size()),lastr(b.size());

	MatrixMultVector<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&z[0]),thrust::raw_pointer_cast(&r[0]),prec.dev_matrix,prec.preA,z.size());
	p = z;
	std::vector<Scalar> tempp{z.begin(),z.end()};
	iter = 0;
	thrust::transform(thrust::device,r.begin(),r.end(),r.begin(),temp.begin(),thrust::multiplies<Scalar>());
	norm = 1000;
	double normb = b.norm1();

	thrust::transform(thrust::device,r.begin(),r.end(),z.begin(),temp.begin(),thrust::multiplies<Scalar>());
	rr0 = thrust::reduce(thrust::device,temp.begin(),temp.end());

	while(iter < limit && norm > tolerance * normb)
	{
		MatrixMultVector<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&w[0]),thrust::raw_pointer_cast(&p[0]),cspm.dev_matrix,cspm.preA,Ap.size());

		std::vector<Scalar> tempp{w.begin(),w.end()};
		thrust::transform(thrust::device,p.begin(),p.end(),w.begin(),temp.begin(),thrust::multiplies<Scalar>());
		alpha = rr0 / thrust::reduce(thrust::device,temp.begin(),temp.end());

		computeX_PCG<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&xx[0]),thrust::raw_pointer_cast(&p[0]),xx.size(),alpha);
		tempp = {xx.begin(),xx.end()};

		lastr = r;
		computeR_PCG<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&r[0]),thrust::raw_pointer_cast(&w[0]),r.size(),alpha);
		tempp = {r.begin(),r.end()};

		MatrixMultVector<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&z[0]),thrust::raw_pointer_cast(&r[0]),prec.dev_matrix,prec.preA,z.size());
		tempp = {z.begin(),z.end()};


		thrust::transform(thrust::device,r.begin(),r.end(),lastr.begin(),temp.begin(),thrust::minus<Scalar>());
		thrust::transform(thrust::device,temp.begin(),temp.end(),z.begin(),temp.begin(),thrust::multiplies<Scalar>());
		rr1 = thrust::reduce(thrust::device,temp.begin(),temp.end());

		beta = rr1 / rr0;
		thrust::transform(thrust::device,r.begin(),r.end(),z.begin(),temp.begin(),thrust::multiplies<Scalar>());
		rr0 = thrust::reduce(thrust::device,temp.begin(),temp.end());
		// rr0 = rr1;

		computeP_PCG<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&p[0]),thrust::raw_pointer_cast(&z[0]),p.size(),beta);
		tempp = {p.begin(),p.end()};

		iter++;
		thrust::transform(thrust::device,r.begin(),r.end(),r.begin(),temp.begin(),thrust::multiplies<Scalar>());
		norm = thrust::reduce(thrust::device,temp.begin(),temp.end());
		norm = std::sqrt(norm);
	}
	x.setvalues({xx.begin(),xx.end()});
}

void BiCGSTAB(const SymetrixSparseMatrix& A,Vector& x,const Vector& b,double tolerance,int limit,int& iter,double& norm)
{
	int bs = (A.get_row() / 32 * 32 == A.get_row()) ?  A.get_row() / 32 : A.get_row() / 32 + 1;
	dim3 blockSize(bs);
	dim3 threadSize(32);
	Scalar rho0,w,alpha,rho1;
	rho0 = w = alpha = 1.0;

	CudaSPMatrix cspm(A.get_row(),A.get_col(),A);
	thrust::device_vector<Scalar> r(b.begin(),b.end()),xx(b.size()),r_hat = r,v(b.size()),p(b.size()),s(b.size()),t(b.size()),temp(b.size());
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
		MatrixMultVector<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&v[0]),thrust::raw_pointer_cast(&p[0]),cspm.dev_matrix,cspm.preA,v.size());

		thrust::transform(thrust::device,r_hat.begin(),r_hat.end(),v.begin(),temp.begin(),thrust::multiplies<Scalar>());
		alpha = rho1 / thrust::reduce(thrust::device,temp.begin(),temp.end());

		computeS<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&s[0]),thrust::raw_pointer_cast(&r[0]),thrust::raw_pointer_cast(&v[0]),s.size(),alpha);
		MatrixMultVector<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&t[0]),thrust::raw_pointer_cast(&s[0]),cspm.dev_matrix,cspm.preA,t.size());


		thrust::transform(thrust::device,s.begin(),s.end(),t.begin(),temp.begin(),thrust::multiplies<Scalar>());
		w = thrust::reduce(thrust::device,temp.begin(),temp.end());
		thrust::transform(thrust::device,t.begin(),t.end(),t.begin(),temp.begin(),thrust::multiplies<Scalar>());
		w = w / thrust::reduce(thrust::device,temp.begin(),temp.end());


		// thrust::transform(thrust::device,r_hat.begin(),r_hat.end(),t.begin(),temp.begin(),thrust::multiplies<Scalar>());
		// rho1 = -w * thrust::reduce(thrust::device,temp.begin(),temp.end());
		computeX<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&xx[0]),thrust::raw_pointer_cast(&p[0]),thrust::raw_pointer_cast(&s[0]),xx.size(),alpha,w);
		computeR<<<blockSize,threadSize>>>(thrust::raw_pointer_cast(&r[0]),thrust::raw_pointer_cast(&s[0]),thrust::raw_pointer_cast(&t[0]),r.size(),w);
		iter++;
		thrust::transform(thrust::device,r.begin(),r.end(),r.begin(),temp.begin(),thrust::multiplies<Scalar>());
		norm = thrust::reduce(thrust::device,temp.begin(),temp.end());
		norm = std::sqrt(norm);
	}
	x.setvalues({xx.begin(),xx.end()});
}