#ifndef GEANTV_VECTORPHYSICSTYPES_H
#define GEANTV_VECTORPHYSICSTYPES_H

#include <GV/VecCore/VecCore>

using PhysVecBackend       = vecCore::backend::VcVector;
using PhysDV               = vecCore::backend::VcVector::Double_v;
using PhysI32V             = vecCore::backend::VcVector::Int32_v;
using PhysI64V             = vecCore::backend::VcVector::Int64_v;
using PhysDM               = vecCore::Mask<PhysDV>;
using PhysDI               = vecCore::Index<PhysDV>;
constexpr int kPhysDVWidth = (int)vecCore::VectorSize<PhysDV>();
constexpr int kPhysDVAlign = (int)vecCore::VectorSize<PhysDV>() * sizeof(double);

#endif // GEANTV_VECTORPHYSICSTYPES_H
