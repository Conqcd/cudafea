//
// Constraint class
//
//  Constraint.cpp
//  
//

#include "Constraint.hpp"


template <typename T> Constraint<T>::Constraint() : cx(1), cy(1), cz(1), node(0), preserve(1) { }

template <typename T> Constraint<T>::Constraint(bool CX, bool CY, bool CZ) : cx(CX), cy(CY), cz(CZ), node(0), preserve(1) { }

template <typename T> Constraint<T>::~Constraint() { }

//=================================================

template class Constraint<xyzType>; 



