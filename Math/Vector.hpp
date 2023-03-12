#pragma once
#include "../Common.hpp"
#include<vector>
class Matrix;
class Vector
{
private:
    std::vector<Scalar> m_Vec;
public:
    Vector(idxType length);
    Vector() = default;
    ~Vector();

    idxType size()const {return m_Vec.size();};

    void reset(idxType);
    void destroy();
    void fill(double);
    Vector& scale(double s);
    void setvalues(const std::vector<idxType>&,const std::vector<Scalar>&);
    void setvalue(idxType,Scalar);
    void setvalues(const std::vector<Scalar>&);
    double norm1()const;

    std::vector<Scalar> generateScalar()const;
    Vector& operator+=(const Vector& v);
    Vector& operator-=(const Vector& v);
    Scalar operator[](unsigned int index)const {return m_Vec[index];}

    inline auto begin()const {return m_Vec.begin();}
    inline auto end()const {return m_Vec.end();}
    friend Vector operator*(const Matrix& m, const Vector& vec);
};

inline Vector operator+(const Vector& v1,const Vector& v2)
{
    Vector v = v1;
    return v += v2;
}

inline Vector operator-(const Vector& v1,const Vector& v2)
{
    Vector v = v1;
    return v -= v2;
}