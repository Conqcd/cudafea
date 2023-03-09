#pragma once

#include "../Math/BiCGSTABSolver.hpp"

struct IndexValue
{
    int colid;
    double value;
};

class CudaVector
{
private:
    int row;
    double* values;
public:
    CudaVector(int _row,const Vector& vec) : row(_row){AllocateData(vec);}
    ~CudaVector();
protected:
    void AllocateData(const Vector&);
};

class CudaSPVector
{
private:
    int row;
    int preA;
    IndexValue* vec;
public:
    CudaSPVector(int _row,const Vector& vec) : row(_row){AllocateData(vec);}
    ~CudaSPVector();
protected:
    void AllocateData(const Vector&);
};

class CudaSPMatrix
{
private:
    int row,col;
    int* preA;
    IndexValue** matrix;
    IndexValue** dev_matrix;
public:
    CudaSPMatrix(int _row,int _col,const SymetrixSparseMatrix& mat) : row(_row),col(_col){AllocateData(mat);}
    ~CudaSPMatrix();
protected:
    void AllocateData(const SymetrixSparseMatrix&);
};

void BiCGSTAB(const SymetrixSparseMatrix& A,Vector& x,const Vector& b);