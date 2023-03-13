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
    double      yx, yy, yz; // strain
    double      fx = 0, fy = 0, fz = 0, fxy = 0, fzy = 0, fzx = 0; // stress
    double von_mises = 0;
    int ct = 0;
    idxType     idx;        // index
    Constraint<T>*    constraint; // constraint pointer
    Element<T>*    elems[NODES_PER_ELEMENT];
    Node();
    ~Node();
};


#endif /* NODE_H */
