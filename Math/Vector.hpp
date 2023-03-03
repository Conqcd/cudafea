#pragma once
#include "../Common.hpp"
#include<vector>
class Vector
{
private:
    idxType m_length;
    std::vector<double> m_Vec;
public:
    Vector(idxType length);
    ~Vector();
};
