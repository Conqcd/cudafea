//
// Material class
//
//
//  Material.h
//
//

#ifndef MATERIAL_H
#define MATERIAL_H

#include "Common.h"
#include "LocalStiffnessMatrix.h"
#include "petscsys.h"

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

