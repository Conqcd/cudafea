#include "Vector.hpp"
#include <cmath>
#include <assert.h>

Vector::Vector(idxType length)
{
    m_Vec.resize(length);
}

Vector::~Vector()
{
}

void Vector::reset(idxType length)
{
    m_Vec.resize(length);
}

void Vector::destroy()
{

}

void Vector::fill(double value)
{
    std::fill(m_Vec.begin(),m_Vec.end(),value);
}

void Vector::setvalues(const std::vector<idxType>& idvec,const std::vector<Scalar>& valvec)
{
    assert(idvec.size() == valvec.size());
    int i = 0;
    for (auto id:idvec)
    {
        m_Vec[id] = valvec[i++];
    }
}

void Vector::setvalues(const std::vector<Scalar>& valvec)
{
    m_Vec = valvec;
}

double Vector::norm1()const
{
    double value = 0;
    for (auto v:m_Vec)
    {
        value += v * v;
    }
    return std::sqrt(value);
    // return 0;
}

std::vector<Scalar> Vector::generateScalar()const
{
    return m_Vec;
}