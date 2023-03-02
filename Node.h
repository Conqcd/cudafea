//
// Node class
//

#ifndef NODE_H
#define NODE_H

#include "Common.h"
#include "petscsys.h"

template <typename T> class Element;
template <typename T> class Constraint;

template<typename T> class Node
{
public:
    T           x,y,z;      // coordinates
    PetscScalar dx, dy, dz; // displacements
    idxType     idx;        // index
    Constraint<T>*    constraint; // constraint pointer
    Element<T>*    elems[NODES_PER_ELEMENT];
    Node();
    ~Node();


};


#endif /* NODE_H */

