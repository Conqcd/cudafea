#include "Matrix.hpp"
#include <cmath>
#include <assert.h>

DenseMatrix::DenseMatrix()
{
}

DenseMatrix::DenseMatrix(idxType row,idxType col)
        : Matrix(row,col)
{
    m_Mat.resize(m_row);
    for(auto& row:m_Mat)
        row.resize(col);
}

DenseMatrix::~DenseMatrix()
{
}

void DenseMatrix::reset(idxType row,idxType col)
{
    m_Mat.clear();
    m_row = row;
    m_col = col;
    m_Mat.resize(m_row);
    for(auto& row:m_Mat)
        row.resize(col);
}

void DenseMatrix::insert(idxType row,idxType col,double value)
{
    m_Mat[row][col] = value;
}

void DenseMatrix::add(idxType row,idxType col,double value)
{
    m_Mat[row][col] += value;
}

void DenseMatrix::scale(double s)
{
    for (auto& row:m_Mat)
        for (auto& col:row)
            col *= s;
}

void DenseMatrix::mult(const Matrix& m)
{
    assert(m_col == m.get_row());
    decltype(m_Mat) mat(m_row);
    
    for (int i = 0; i < m_row; i++)
    {
        mat[i].resize(m.get_col());
        for (int j = 0; j < mat[i].size(); j++)
        {
            for (int k = 0; k < m_col; k++)
            {
                mat[i][j] += m_Mat[i][k] * m.index(k,j);
            }
        }
    }
    m_col = m.get_col();
    m_Mat = mat;
}

void DenseMatrix::mult(const Matrix& m1,const Matrix& m2)
{
    assert(m1.get_col() == m2.get_row());
    m_row = m1.get_row();
    m_col = m2.get_col();
    decltype(m_Mat) mat(m1.get_row());
    for (int i = 0; i < m_row; i++)
    {
        mat[i].resize(m2.get_col());
        for (int j = 0; j < mat[i].size(); j++)
        {
            for (int k = 0; k < m1.get_col(); k++)
            {
                mat[i][j] += m1.index(i,k) * m2.index(k,j);
            }
        }
    }
    m_Mat = mat;
}

void DenseMatrix::AXPY(double a,const Matrix& x)
{
    assert(m_row == x.get_row() && m_col == x.get_col());
    int c = 0,r = 0;
    for(auto& row:m_Mat)
    {
        c = 0;
        for(auto&col : row)
        {
            col += a * x.index(r,c);
            c++;
        }
        r++;
    }
}

void DenseMatrix::destroy()
{
        
}

void DenseMatrix::insertValues(const std::vector<idxType>& rowid,const std::vector<idxType>& colid,const std::vector<Scalar>& values)
{
    int id = 0;
    for(auto& row:rowid)
        for(auto& col:colid)
        {
            assert(row < m_row && row >= 0 && col < m_col && col >= 0);
            m_Mat[row][col] = values[id++];
        }
}

void DenseMatrix::PreAllocation(idxType num)
{
    
}

double DenseMatrix::index(idxType row,idxType col)const
{
    return m_Mat[row][col];
}



SymetrixSparseMatrix::SymetrixSparseMatrix()
                    : preA(0)
{

}

SymetrixSparseMatrix::SymetrixSparseMatrix(idxType row,idxType col)
        : Matrix(row,col) ,preA(0)
{
    m_Mat.resize(m_row);
}

SymetrixSparseMatrix::~SymetrixSparseMatrix()
{
}

void SymetrixSparseMatrix::reset(idxType row,idxType col)
{
    m_Mat.clear();
    m_row = row;
    m_col = col;
    m_Mat.resize(m_row);
}

void SymetrixSparseMatrix::insert(idxType row,idxType col,double value)
{
    if(m_Mat[row].count(col) == 0)
    {
        assert(m_Mat[row].size() < preA);
    }
    m_Mat[row][col] = value;
}

void SymetrixSparseMatrix::add(idxType row,idxType col,double value)
{
    if(m_Mat[row].count(col) == 0)
    {
        assert(m_Mat[row].size() < preA);
        m_Mat[row][col] = value;
    }else
    {
        m_Mat[row][col] += value;
    }
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
    assert(m_col == m.get_row());
    // decltype(m_Mat) mat(m_row);

    // for (int i = 0; i < m_row; i++)
    // {
    //     mat.resize(m.get_col());
    //     for (int j = 0; j < mat[i].size(); j++)
    //     {
    //         for (int k = 0; k < mat[i].size(); k++)
    //         {
    //             mat[i][j] += m_Mat[i][k] * m.index(k,j);
    //         }
    //     }
    // }
    // m_Mat = mat;
}

void SymetrixSparseMatrix::mult(const Matrix& m1,const Matrix& m2)
{
    assert(m1.get_col() == m2.get_row());
}

void SymetrixSparseMatrix::AXPY(double a,const Matrix& x)
{
    assert(m_row == x.get_row() && m_col == x.get_col());
}

void SymetrixSparseMatrix::destroy()
{

}

void SymetrixSparseMatrix::insertValues(const std::vector<idxType>& rowid,const std::vector<idxType>& colid,const std::vector<Scalar>& values)
{
    int id = 0;
    for(auto& row:rowid)
        for(auto& col:colid)
        {
            assert(row < m_row && row >= 0 && col < m_col && col >= 0);
            bool exist = false;
            if(m_Mat[row].count(col) == 0)
            {
                assert(m_Mat[row].size() <= preA);
            }
            if(values[id] != 0)
                m_Mat[row][col] = values[id++];
            else
                id++;
    }
}

void SymetrixSparseMatrix::PreAllocation(idxType num)
{
    preA = std::min(num,m_col);
    // for (auto& row:m_Mat)
    // {
    //     row.resize(preA);
    // }
}

double SymetrixSparseMatrix::index(idxType row,idxType col)const
{
    if(m_Mat[row].count(col) == 0)
        return 0;
    return m_Mat[row].at(col);
}
SymetrixSparseMatrix SymetrixSparseMatrix::inverse()const
{
    SymetrixSparseMatrix mat(m_row,m_col);
    mat.PreAllocation(preA);
    for (int i = 0; i < m_row; i++)
    {
        mat.insert(i,i,1 / index(i,i));
        for (int j = 0; j < i; j++)
        {
            Scalar s = 0;
            for (int k = j; k < i; k++)
            {
                s = s + index(i,k) * mat.index(k,j);
            }
            mat.insert(i,j,-s * mat.index(i, i));
        }
    }
    return mat;
}

SymetrixSparseMatrix SymetrixSparseMatrix::ichol()const
{
    auto mat = *this;
    
    mat.PreAllocation(preA);
    for (int k = 0; k < m_row; k++)
    {
        mat.insert(k,k,std::sqrt(mat.index(k,k)));
        
        for (int i = k + 1; i < m_col; i++)
            if(mat.index(i,k) != 0)
                mat.insert(i,k,mat.index(i,k) / mat.index(k,k));
        for (int j = k + 1; j < m_col; j++)
            for (int i = j; i < m_col; i++)
                if (mat.index(i,j) != 0)
                    mat.insert(i,j,mat.index(i,j) - mat.index(i,k) * mat.index(j,k));
    }
    for (int i = 0; i < m_row; i++)
        for (int j = i + 1; j < m_col; j++)
            mat.insert(i,j,0);
    return mat;
}
    
Vector operator*(const Matrix& matrix,const Vector& vec)
{
    Vector vec2(matrix.get_row());
    for (int i = 0; i < matrix.get_row(); i++)
    {
        double v = 0;
        for (int j = 0; j < matrix.get_col(); j++)
        {
            v += matrix.index(i,j) * vec[j];
        }
        
        vec2.setvalue(i,v);
    }
    
    return vec2;
}

namespace Math
{

DenseMatrix transpose(const DenseMatrix& m)
{
    DenseMatrix res(m.get_col(),m.get_row());
    for (int i = 0; i < res.get_row(); i++)
    {
        for (int j = 0; j < res.get_col(); j++)
        {
            res.insert(i,j,m.index(j,i));
        }
    }
    return res;
}


} // namespace Math
