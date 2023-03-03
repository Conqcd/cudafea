#include "Matrix.hpp"

Matrix::Matrix()
        : m_row(0) , m_col(0)
{
}

Matrix::Matrix(idxType row,idxType col)
        : m_row(row) , m_col(col)
{
    m_Mat.resize(m_row);
}

Matrix::~Matrix()
{
}

void Matrix::reset(idxType row,idxType col)
{

}

void Matrix::insert(idxType row,idxType col,double value)
{

}

void Matrix::add(idxType row,idxType col,double value)
{

}

void Matrix::scale(double s)
{

}

void Matrix::mult(const Matrix& m)
{

}

void Matrix::mult(const Matrix& m1,const Matrix& m2)
{

}

void Matrix::AXPY(double a,const Matrix& x)
{

}

void Matrix::destroy()
{
        
}

namespace Math
{
Matrix transpose(const Matrix& m)
{
    return {};
}
    
} // namespace Math
