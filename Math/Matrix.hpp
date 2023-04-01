#pragma once
#include "../Common.hpp"
#include "Vector.hpp"
#include <assert.h>
#include<vector>
#include<map>

template<typename T>
using IdVal = std::pair<idxType,T>;


class Vector;
class Matrix
{
protected:
    idxType m_row,m_col;
public:
    Matrix() : m_row(0),m_col(0){}
    Matrix(idxType row,idxType col):m_row(row),m_col(col){}
    virtual~Matrix(){}

    virtual void reset(idxType row,idxType col) = 0;
    virtual void insert(idxType row,idxType col,double value) = 0;
    virtual void add(idxType row,idxType col,double value) = 0;
    virtual void scale(double s) = 0;
    virtual void mult(const Matrix& m) = 0;
    virtual void mult(const Matrix& m1,const Matrix& m2) = 0;
    virtual void AXPY(double,const Matrix&) = 0;
    virtual void destroy() = 0;
    virtual void insertValues(const std::vector<idxType>&,const std::vector<idxType>&,const std::vector<Scalar>&) = 0;
    virtual void PreAllocation(idxType) = 0;
    virtual Scalar index(idxType,idxType)const = 0;
    virtual Scalar index(idxType,idxType) = 0;


    inline idxType get_row()const {return m_row;}
    inline idxType get_col()const {return m_col;}
    virtual double operator[](unsigned int index)const = 0;
    friend Vector operator*(const Matrix& m, const Vector& vec);
};
class DenseMatrix :public Matrix
{
private:
    std::vector<std::vector<double>> m_Mat;
public:
    DenseMatrix();
    DenseMatrix(idxType row,idxType col);
    virtual~DenseMatrix();

    virtual void reset(idxType row,idxType col) override;
    virtual void insert(idxType row,idxType col,double value) override;
    virtual void add(idxType row,idxType col,double value) override;
    virtual void scale(double s) override;
    virtual void mult(const Matrix& m) override;
    virtual void mult(const Matrix& m1,const Matrix& m2) override;
    virtual void AXPY(double,const Matrix&) override;
    virtual void destroy() override;
    virtual void insertValues(const std::vector<idxType>&,const std::vector<idxType>&,const std::vector<Scalar>&) override;
    virtual void PreAllocation(idxType) override;
    virtual Scalar index(idxType row,idxType col)const override;
    virtual Scalar index(idxType row,idxType col)override;

    virtual Scalar operator[](unsigned int index)const override
    {
        return m_Mat[index / m_col][index % m_col];
    }
};

class SymetrixSparseMatrix : public Matrix
{
private:
    std::vector<std::vector<std::pair<idxType,Scalar>>> m_Mat;
    idxType preA;
    idxType count;
public:
    SymetrixSparseMatrix();
    SymetrixSparseMatrix(idxType row,idxType col);
    virtual~SymetrixSparseMatrix();

    virtual void reset(idxType row,idxType col)override;
    virtual void insert(idxType row,idxType col,Scalar value)override;
    virtual void add(idxType row,idxType col,Scalar value)override;
    virtual void scale(Scalar s)override;
    virtual void mult(const Matrix& m)override;
    virtual void mult(const Matrix& m1,const Matrix& m2)override;
    virtual void AXPY(Scalar,const Matrix&)override;
    virtual void destroy()override;
    virtual void insertValues(const std::vector<idxType>&,const std::vector<idxType>&,const std::vector<Scalar>&)override;
    virtual void PreAllocation(idxType)override;
    virtual Scalar index(idxType row,idxType col)const override;
    virtual Scalar index(idxType row,idxType col)override;

    SymetrixSparseMatrix inverse_lowertri()const;
    SymetrixSparseMatrix transpose()const;
    SymetrixSparseMatrix ichol()const;
    SymetrixSparseMatrix SSORAI()const;

    void SolveTriL(Vector& x,const Vector& b);
    void SolveTriU(Vector& x,const Vector& b);

    inline const auto& getRow(idxType row)const {
        assert(row < m_row && row >= 0);
        return m_Mat[row];
    }
    inline auto getCol(idxType row,idxType col)const{
        std::pair<idxType,Scalar> p;
        p.first = col;
        return std::lower_bound(m_Mat[row].begin(),m_Mat[row].end(),p,[](const auto& kv1,const auto& kv2) -> bool {return kv1.first < kv2.first;});
    }
    inline auto getCol(idxType row,idxType col){
        std::pair<idxType,Scalar> p;
        p.first = col;
        return std::lower_bound(m_Mat[row].begin(),m_Mat[row].end(),p,[](const auto& kv1,const auto& kv2) -> bool {return kv1.first < kv2.first;});
    }
    virtual Scalar operator[](unsigned int id)const override
    {
        return index(id / m_col,id % m_col);
    }
    // SymetrixSparseMatrix operator *(const Matrix& m);
};

Vector operator*(const Matrix& matrix,const Vector& vec);
namespace Math
{

DenseMatrix transpose(const DenseMatrix& m);


}