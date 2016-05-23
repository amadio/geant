#ifndef RANDOM_KERNEL_H
#define RANDOM_KERNEL_H

#include "base/VPGlobal.h"

namespace vecphys {

bool GUCurand_Init(Random_t* randomStates,
		   unsigned long seed,
		   int blocksPerGrid,
		   int threadsPerBlock);

} // end namespace vecphys

#endif
