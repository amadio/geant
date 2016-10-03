// Toy physics processes for our propagator prototype. Currently including:
// - single scattering as a discrete process
// - energy loss as continuous process
// - generic interaction as discrete process, producing secondaries

#include "PhysicsProcessOld.h"

// pin virtual table.
PhysicsProcessOld::~PhysicsProcessOld() {}
