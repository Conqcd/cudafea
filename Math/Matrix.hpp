#pragma once
#include "../Common.hpp"
#include<vector>

template<typename T>
using IdVal = std::pair<idxType,T>;

//defaut symatric
class Matrix
{
private:
    idxType m_row,m_col;
    std::vector<std::vector<IdVal<double>>> m_Mat;
public:
    Matrix();
    Matrix(idxType row,idxType col);
    ~Matrix();

    void reset(idxType row,idxType col);
    void insert(idxType row,idxType col,double value);
    void add(idxType row,idxType col,double value);
    void scale(double s);
    void mult(const Matrix& m);
    void mult(const Matrix& m1,const Matrix& m2);
    void AXPY(double,const Matrix&);
    void destroy();
    void insertValues(const std::vector<idxType>&,const std::vector<idxType>&,const std::vector<Scalar>&);
    void PreAllocation(idxType);


    idxType get_row()const {return m_row;}
    idxType get_col()const {return m_col;}
    double operator[](unsigned int index)const
    {
        return 0;
    }
};

namespace Math
{
Matrix transpose(const Matrix& m);

}