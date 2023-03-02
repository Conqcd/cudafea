//
// Element class
//
//  Element.h
//  
//

#ifndef ELEMENT_H
#define ELEMENT_H

#include "Common.h"

template <typename T> class Node;

template <typename T>
class Element{
public :
    T         x, y, z;
    idxType   idx;
    midxType  material; // material idx 
    Node<T> * nodes[NODES_PER_ELEMENT];
    Element();
    Element(midxType M);
    ~Element();
};

#endif /* ELEMENT_H */