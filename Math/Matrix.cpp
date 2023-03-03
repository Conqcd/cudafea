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