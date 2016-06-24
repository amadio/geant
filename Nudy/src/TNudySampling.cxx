// 	Selection of secondary particle energy, angle 
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: June 22, 2016

#include <iostream>
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
TNudySampling::TNudySampling(){}
//______________________________________________________________________________
TNudySampling::TNudySampling(const char *rENDF)
{
  ///*
  recoPoint = new TNudyEndfRecoPoint();
  recoAng = new TNudyEndfAng ();
  recoPoint->SetsigPrecision(1E-3);
  recoPoint->GetData(rENDF);
  double sigmaTotal = recoPoint->GetSigmaTotal(20.0);
  std::cout <<" sigma Total "<< sigmaTotal << std::endl;
  /*
  for (unsigned int i = 0; i < recoPoint->sigmaUniOfMts.size(); i++) {
    double sigmaPartial = recoPoint->GetSigmaPartial(i, 20.0);
    std::cout <<"MT = "<< recoPoint->MtValues[0][i] <<" sigma Partial = "<< sigmaPartial << std::endl;
     std::cout<<" cos "<< recoPoint->GetCos4(recoPoint->MtValues[0][i], 20.0) << std::endl;
   std::cout<<" energy "<< recoPoint->GetEnergy5(recoPoint->MtValues[0][i], 20.0) << std::endl;
  }
  */
  //*/
} 
TNudySampling::~TNudySampling(){}
