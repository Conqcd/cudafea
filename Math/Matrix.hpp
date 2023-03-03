#pragma once
#include "../Common.hpp"
#include<vector>
class Matrix
{
private:
    idxType m_row,m_col;
    std::vector<std::vector<double>> m_Mat;
public:
    Matrix();
    Matrix(idxType row,idxType col);
    ~Matrix();

    void reset(idxType row,idxType col);
    void insert(idxType row,idxType col,double value);
    void scale(double s);
    void mult(const Matrix& m);
};
