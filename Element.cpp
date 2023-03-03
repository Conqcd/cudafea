//
// Element class
//
// Element.cpp
//  
//

#include "Element.hpp"

//====================================================

template <typename T> Element<T>::Element()
{
    for(int kk = 0; kk < NODES_PER_ELEMENT; ++kk) nodes[kk] = 0;
}

template <typename T> Element<T>::Element(midxType M) : material(M)
{
    for(int kk = 0; kk < NODES_PER_ELEMENT; ++kk) nodes[kk] = 0;
}

template <typename T> Element<T>::~Element() { }

//====================================================

template class Element<xyzType>;

