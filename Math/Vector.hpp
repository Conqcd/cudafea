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
    double norm1()const;

    std::vector<Scalar> generateScalar()const;
    Scalar operator[](unsigned int index)const {return m_Vec[index];}

    inline auto begin()const {return m_Vec.begin();}
    inline auto end()const {return m_Vec.end();}
};
