#include "Vector.hpp"

Vector::Vector(idxType length)
        : m_length(length)
{
    m_Vec.resize(m_length);
}

Vector::~Vector()
{
}
