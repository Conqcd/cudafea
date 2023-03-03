#pragma once
#include "../Common.hpp"
#include<vector>
class Vector
{
private:
    std::vector<double> m_Vec;
public:
    Vector(idxType length);
    Vector() = default;
    ~Vector();

    idxType size()const {return m_Vec.size();};

    void reset(idxType);
    void destroy();
};
