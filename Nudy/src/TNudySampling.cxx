// 	Selection of secondary particle energy, angle 
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: June 22, 2016

#include "TNudyEndfNuPh.h"
#include "TNudyEndfDoppler.h"
#include "TNudyEndfAng.h"
#include "TNudyEndfEnergy.h"
#include "TNudyEndfEnergyAng.h"
#include "TNudyEndfFissionYield.h"
#include "TNudyCore.h"
#include "TNudyEndfRecoPoint.h"
#include "TNudySampling.h"

#ifdef USE_ROOT
ClassImp(TNudySampling)
#endif
//______________________________________________________________________________
TNudySampling::TNudySampling()
{
//double sigmat = recoPoint.GetSigmaTotal(1.0);
//std::cout << sigmat << std::endl;

  
}

TNudySampling::~TNudySampling(){}
