#include "Vector.hpp"

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

}

void Vector::setvalues(const std::vector<idxType>& idvec,const std::vector<Scalar>& valvec)
{

}

std::vector<Scalar> Vector::generateScalar()
{
    return m_Vec;
}