// 	Selection of secondary particle energy, angle 
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: June 22, 2016
#include <iostream>
//#include "TNudyEndfDoppler.h"
//#include "TNudyEndfAng.h"
//#include "TNudyEndfEnergy.h"
//#include "TNudyEndfEnergyAng.h"
//#include "TNudyEndfFissionYield.h"
//#include "TNudyCore.h"
#include "TNudyEndfRecoPoint.h"
#include "TNudySampling.h"
#ifdef USE_ROOT
#include "TRandom.h"
#endif

#ifdef USE_ROOT
ClassImp(TNudySampling)
#endif

TNudySampling::TNudySampling(){}
//------------------------------------------------------------------------------------------------------
TNudySampling::TNudySampling(TNudyEndfRecoPoint *recoPoint){
  fRnd = new TRandom();
  std::cout <<"sampling sigmaTotal "<< recoPoint->GetSigmaTotal(0,20) << std::endl;
  std::cout <<"sampling sigmaPartial total "<< recoPoint->GetSigmaPartial(0,0,20) << std::endl;
  std::cout <<"sampling sigmaPartial elstic "<< recoPoint->GetSigmaPartial(0,1,20) << std::endl;
  // determining reaction type from element;
  int elemId = 0;
  //double density = 1;
  //double mass = 1;
  //double charge = 1;
  //double avg = 6.022E23;
  //double ro = avg * density / mass;
  double kineticE = 10; 
  int isel = 0;
  std::vector<double> crs;
  int enemax = recoPoint->eneUni[elemId].size(); 
  int counter =0;
  do
  {
    double sum1 = 0;
    kineticE = fRnd->Uniform(1) * recoPoint->eneUni[elemId][enemax-1];
    for(unsigned int crsp = 0; crsp < recoPoint->MtValues[elemId].size(); crsp++){
      crs.push_back(recoPoint->GetSigmaPartial(elemId,crsp,kineticE)/recoPoint->GetSigmaTotal(elemId,kineticE));
    }  
    double rnd1 = fRnd->Uniform(1);    
    for(unsigned int crsp = 0; crsp < recoPoint->MtValues[elemId].size(); crsp++){
      sum1 += crs[crsp];
      if(rnd1 <= sum1){
	isel = crsp;
	break;
      }
    }
    //std::cout <<"selected reaction MT  "<< recoPoint->MtValues[elemId][isel] <<"  "<< kineticE << std::endl;
    int MT = recoPoint->MtValues[elemId][isel];
    //double cosT = recoPoint->GetCos4(elemId, MT, kineticE);
    //std::cout <<"counter "<<counter <<" cos "<< cosT << std::endl;
    double secEnergy = recoPoint->GetEnergy5(elemId, MT, kineticE);
    //std::cout <<"counter "<<counter <<" energy "<< secEnergy << std::endl;
    crs.clear();
    counter++;
  }while(counter < 100);
  std::cout<<" end sampling "<< std::endl;
}
//------------------------------------------------------------------------------------------------------

TNudySampling::~TNudySampling(){}
