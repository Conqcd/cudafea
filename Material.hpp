//
// Material class
//
//
//  Material.h
//
//

#ifndef MATERIAL_H
#define MATERIAL_H

#include "Common.hpp"
#include "LocalStiffnessMatrix.hpp"
// #include "petscsys.h"

 class Material
{
    
public:
    midxType        idx;
    double          youngsm;
    double          poissonsr;
    bool            remodel;
    LocalStiffnessMatrix lsm;
    
    Material();
    Material(midxType id, double ym, double pr, LocalStiffnessMatrix mat, bool rf);
    ~Material();
};


#endif /* MATERIAL_H */

