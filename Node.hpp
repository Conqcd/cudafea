//
// Node class
//

#ifndef NODE_H
#define NODE_H

#include "Common.hpp"

template <typename T> class Element;
template <typename T> class Constraint;

template<typename T> class Node
{
public:
    T           x,y,z;      // coordinates
    double      dx, dy, dz; // displacements
    idxType     idx;        // index
    Constraint<T>*    constraint; // constraint pointer
    Element<T>*    elems[NODES_PER_ELEMENT];
    Node();
    ~Node();
};


#endif /* NODE_H */
