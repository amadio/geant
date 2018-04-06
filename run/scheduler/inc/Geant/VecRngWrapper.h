#ifndef GEANTV_VECRNGWRAPPER_H
#define GEANTV_VECRNGWRAPPER_H

#include <GV/VecRng/MRG32k3a.h>

#include "Geant/VectorPhysicsTypes.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

class VecRngWrapper {
public:
  VecRngWrapper()
  {
    void *buff     = vecCore::AlignedAlloc(64, sizeof(vecRng::MRG32k3a<vecCore::backend::Scalar>));
    mrg32k3aScalar = new (buff) vecRng::MRG32k3a<vecCore::backend::Scalar>;

    buff        = vecCore::AlignedAlloc(64, sizeof(vecRng::MRG32k3a<PhysVecBackend>));
    mrg32k3aVec = new (buff) vecRng::MRG32k3a<PhysVecBackend>;

    mrg32k3aScalar->Initialize();
    mrg32k3aVec->Initialize();
  }
  ~VecRngWrapper()
  {
    // mrg32k3aScalar->vecRng::~MRG32k3a<vecCore::backend::Scalar>();
    vecCore::AlignedFree(mrg32k3aScalar);
    // mrg32k3aVec->vecRng::~MRG32k3a();
    vecCore::AlignedFree(mrg32k3aVec);
  }
  void uniform_array(size_t n, double *array, const double min = 0., const double max = 1.)
  {
    for (size_t i = 0; i < n; ++i) {
      array[i] = uniform(min, max);
    }
  }

  double uniform() { return mrg32k3aScalar->Uniform<vecCore::backend::Scalar>(); }
  double uniform(double a, double b) { return a + (b - a) * uniform(); }
  PhysDV uniformV() { return mrg32k3aVec->Uniform<PhysVecBackend>(); }

  double Gauss(double mean, double sigma) { return mrg32k3aScalar->Gauss<vecCore::backend::Scalar>(mean, sigma); }

private:
  vecRng::MRG32k3a<vecCore::backend::Scalar> *mrg32k3aScalar;
  vecRng::MRG32k3a<PhysVecBackend> *mrg32k3aVec;
};
}
}

#endif // GEANTV_VECRNGWRAPPER_H
