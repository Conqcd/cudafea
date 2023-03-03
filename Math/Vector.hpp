#pragma once
#include "../Common.hpp"
#include<vector>
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
    void setvalues(const std::vector<idxType>&,const std::vector<Scalar>&);

    std::vector<Scalar> generateScalar();
};
