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
    std::fill(m_Vec.begin(),m_Vec.end(),value);
}

void Vector::setvalues(const std::vector<idxType>& idvec,const std::vector<Scalar>& valvec)
{
    assert(idvec.size() == valvec.size());
    int i = 0;
    for (auto& id:idvec)
    {
        m_Vec[id] = valvec[id++];
    }
}

std::vector<Scalar> Vector::generateScalar()
{
    return m_Vec;
}