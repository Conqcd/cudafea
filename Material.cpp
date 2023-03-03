//
// Material class
//
//  Material.cpp
//  
//

#include "Material.hpp"

Material::Material() : idx(0), youngsm(0.0), poissonsr(-1.0), remodel(0) { }

Material::Material(midxType id, double ym, double pr, LocalStiffnessMatrix mat, bool rf) : idx(id), youngsm(ym), poissonsr(pr), lsm(mat), remodel(rf) { }

Material::~Material() { }

