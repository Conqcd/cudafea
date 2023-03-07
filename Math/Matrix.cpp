#include "Matrix.hpp"

DenseMatrix::DenseMatrix()
{
}

DenseMatrix::DenseMatrix(idxType row,idxType col)
        : Matrix(row,col)
{
}

DenseMatrix::~DenseMatrix()
{
}

void DenseMatrix::reset(idxType row,idxType col)
{
    m_Mat.resize(m_row);
    m_col = col;
}

void DenseMatrix::insert(idxType row,idxType col,double value)
{

}

void DenseMatrix::add(idxType row,idxType col,double value)
{

}

void DenseMatrix::scale(double s)
{

}

void DenseMatrix::mult(const Matrix& m)
{

}

void DenseMatrix::mult(const Matrix& m1,const Matrix& m2)
{

}

void DenseMatrix::AXPY(double a,const Matrix& x)
{

}

void DenseMatrix::destroy()
{
        
}

void DenseMatrix::insertValues(const std::vector<idxType>& rowid,const std::vector<idxType>& colid,const std::vector<Scalar>& values)
{

}

void DenseMatrix::PreAllocation(idxType num)
{
    for (auto& row:m_Mat)
    {
        row.resize(num);
    }
}

SymetrixSparseMatrix::SymetrixSparseMatrix()
{

}

SymetrixSparseMatrix::SymetrixSparseMatrix(idxType row,idxType col)
        : Matrix(row,col)
{
    m_Mat.resize(m_row);
}

SymetrixSparseMatrix::~SymetrixSparseMatrix()
{
}

void SymetrixSparseMatrix::reset(idxType row,idxType col)
{
    m_Mat.clear();
    m_Mat.resize(m_row);
    m_col = col;
}

void SymetrixSparseMatrix::insert(idxType row,idxType col,double value)
{

}

void SymetrixSparseMatrix::add(idxType row,idxType col,double value)
{

}

void SymetrixSparseMatrix::scale(double s)
{
    for (auto& row:m_Mat)
    {
        for (auto& item:row)
        {
            item.second *= s;
        }
    }
}

void SymetrixSparseMatrix::mult(const Matrix& m)
{

}

void SymetrixSparseMatrix::mult(const Matrix& m1,const Matrix& m2)
{

}

void SymetrixSparseMatrix::AXPY(double a,const Matrix& x)
{

}

void SymetrixSparseMatrix::destroy()
{
        
}

void SymetrixSparseMatrix::insertValues(const std::vector<idxType>& rowid,const std::vector<idxType>& colid,const std::vector<Scalar>& values)
{

}

void SymetrixSparseMatrix::PreAllocation(idxType num)
{
    for (auto& row:m_Mat)
    {
        row.resize(num);
    }
}
namespace Math
{

DenseMatrix transpose(const DenseMatrix& m)
{
    return {};
}
    
} // namespace Math
