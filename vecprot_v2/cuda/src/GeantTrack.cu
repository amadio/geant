
#include "GeantCudaUtils.h"
#include "backend/cuda/Interface.h"

#include "GeantTrack.h"
#include "GeantTaskData.h"
#include "CoprocessorBrokerKernel.h"

void CoprocessorBrokerInitConstant() {}

namespace Geant {
inline namespace cuda {

   template void MakeInstanceAt(GeantTrack_v *addr, unsigned int, int);

} // cuda
} // Geant

namespace vecgeom {
namespace cxx {
   template size_t DevicePtr<Geant::cuda::GeantTaskData>::SizeOf();
   template void DevicePtr<Geant::cuda::GeantTaskData>::ConstructArray(unsigned long, unsigned int, int, unsigned int) const;
   template size_t DevicePtr<Geant::cuda::GeantTrack_v>::SizeOf();
} // cxx
} // vecgeom

