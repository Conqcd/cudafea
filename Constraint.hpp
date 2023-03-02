//
// Constraint class
//  

#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include "Common.h"

template <typename T> class Node;

template <typename T> class Constraint
{
    
public:
    bool      cx, cy, cz;  // 0 if constrained, 1 if not
    double    vx, vy, vz;  // value of constraint
    bool      preserve; // preserve the node 
    // double  dx, dy, dz;  // RHS adjustment due to constraint
    Node<T> * node;        // points to constrained node
    Constraint();
    Constraint(bool CX, bool CY, bool CZ);
    ~Constraint();
};


#endif /* CONSTRAINT_H*/