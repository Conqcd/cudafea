#pragma once
#include "../Common.hpp"
#include<vector>

template<typename T>
using IdVal = std::pair<idxType,T>;

//defaut symatric
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


    idxType get_row()const {return m_row;}
    idxType get_col()const {return m_col;}
    virtual double operator[](unsigned int index)const = 0;
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

    virtual double operator[](unsigned int index)const override
    {
        return 0;
    }

};

class SymetrixSparseMatrix : public Matrix
{
private:
    std::vector<std::vector<IdVal<double>>> m_Mat;
public:
    SymetrixSparseMatrix();
    SymetrixSparseMatrix(idxType row,idxType col);
    virtual~SymetrixSparseMatrix();

    virtual void reset(idxType row,idxType col)override;
    virtual void insert(idxType row,idxType col,double value)override;
    virtual void add(idxType row,idxType col,double value)override;
    virtual void scale(double s)override;
    virtual void mult(const Matrix& m)override;
    virtual void mult(const Matrix& m1,const Matrix& m2)override;
    virtual void AXPY(double,const Matrix&)override;
    virtual void destroy()override;
    virtual void insertValues(const std::vector<idxType>&,const std::vector<idxType>&,const std::vector<Scalar>&)override;
    virtual void PreAllocation(idxType)override;


    virtual double operator[](unsigned int index)const override
    {
        return 0;
    }
};
namespace Math
{

DenseMatrix transpose(const DenseMatrix& m);

}