#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include<thrust/host_vector.h>
#include<thrust/device_vector.h>
#include<thrust/copy.h>
#include<thrust/sort.h>
#include<thrust/execution_policy.h>

#include "../Math/BiCGSTABSolver.hpp"

struct IndexValue
{
    int colid;
    double value;
    // __host__ __device__  IndexValue(){}
    // __host__ __device__  ~IndexValue(){}
    // __host__ __device__  IndexValue(const IndexValue& v){
    //     colid = v.colid;
    //     value = v.colid;
    // }
    // __host__ __device__ IndexValue operator=(const IndexValue& v)
    // {
    //     colid = v.colid;
    //     value = v.colid;
    //     return *this;
    // }
};
class CudaVector
{
private:
    thrust::device_vector<Scalar> values;
public:
    CudaVector(const Vector& vec) {AllocateData(vec);}
    ~CudaVector();
protected:
    void AllocateData(const Vector&);
};

class CudaSPVector
{
private:
    thrust::device_vector<idxType> row;
    thrust::device_vector<Scalar> values;
public:
    CudaSPVector(const Vector& vec) {AllocateData(vec);}
    ~CudaSPVector();
protected:
    void AllocateData(const Vector&);
};

class CudaSPMatrix
{
private:
    int row,col;
    IndexValue** matrix;
public:
    int* preA;
    IndexValue** dev_matrix;
    CudaSPMatrix(int _row,int _col,const SymetrixSparseMatrix& mat) : row(_row),col(_col){AllocateData(mat);}
    ~CudaSPMatrix();
protected:
    void AllocateData(const SymetrixSparseMatrix&);
};

void BiCGSTAB(const SymetrixSparseMatrix& A,Vector& x,const Vector& b,double tolerance,int limit,int& iter,double& norm);