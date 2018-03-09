#ifndef GEANTV_VECTORPHYSICSTYPES_H
#define GEANTV_VECTORPHYSICSTYPES_H

#include <GV/VecCore/VecCore>

using PhysVecBackend       = vecCore::backend::VcVector;
using PhysDV               = vecCore::backend::VcVector::Double_v;
using PhysDM               = vecCore::Mask<PhysDV>;
constexpr int kPhysDVWidth = (int)vecCore::VectorSize<PhysDV>();

#endif // GEANTV_VECTORPHYSICSTYPES_H
