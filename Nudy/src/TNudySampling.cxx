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
#include "TNudyCore.h"
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
TNudySampling::TNudySampling(Particle* particle, TNudyEndfRecoPoint *recoPoint){
  fRnd = new TRandom();
  std::cout <<"sampling sigmaTotal "<< recoPoint->GetSigmaTotal(0,20) << std::endl;
  std::cout <<"sampling sigmaPartial total "<< recoPoint->GetSigmaPartial(0,0,20) << std::endl;
  std::cout <<"sampling sigmaPartial elstic "<< recoPoint->GetSigmaPartial(0,1,20) << std::endl;
  // determining reaction type from element;
  int elemId = 0;
  //double density = 1;
  //double charge = 1;
  //double avg = 6.022E23;
  //double ro = avg * density / mass;
  double kineticE = 10; 
  int isel = 0;
  std::vector<double> crs;
  int enemax = recoPoint->eneUni[elemId].size(); 
  int counter =0;
  double cosCM=0, cosLab=0, secEnergyCM=0, secEnergyLab=0;
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
    int MT = recoPoint->MtValues[elemId][isel];
    int MF = recoPoint->GetMt456 ( elemId, MT ) ;
    std::cout <<" MF "<< MF << std::endl;
    // selection of the data from file 4 5 and 6 for angle and energy calculations
    switch (MF) {
      case 4:
	if(MT == 2){
	  cosCM = recoPoint->GetCos4(elemId, MT, kineticE);
	  cosLab = TNudyCore::Instance()->cmToLabElasticCosT(cosCM, particle[elemId].mass);
	  secEnergyLab = TNudyCore::Instance()->cmToLabElasticE(kineticE, cosCM, particle[elemId].mass);
	} else if(MT >= 50 && MT <= 91){
          
	  cosCM = recoPoint->GetCos4(elemId, MT, kineticE);
	  secEnergyCM = recoPoint->GetEnergy5(elemId, MT, kineticE);
	  secEnergyLab = TNudyCore::Instance()->cmToLabInelasticE(secEnergyCM, kineticE, cosCM, particle[elemId].mass);
	  cosLab = TNudyCore::Instance()->cmToLabInelasticCosT(secEnergyLab, secEnergyCM, kineticE, cosCM, particle[elemId].mass);	  
	}
	break;
      case 5:
	secEnergyCM = recoPoint->GetEnergy5(elemId, MT, kineticE);
	secEnergyLab = TNudyCore::Instance()->cmToLabInelasticE(secEnergyCM, kineticE, cosCM, particle[elemId].mass);
	cosLab = TNudyCore::Instance()->cmToLabInelasticCosT(secEnergyLab, secEnergyCM, kineticE, cosCM, particle[elemId].mass);	  
	break;
      case 6:
	switch (recoPoint->GetLaw6(elemId, MT)){
	  case 1:
	    cosLab = recoPoint->GetCos6(elemId, MT, kineticE);
	    secEnergyLab = recoPoint->GetEnergy6(elemId, MT, kineticE);
	    break;
	  case 2:
	    cosCM = recoPoint->GetCos6(elemId, MT, kineticE);
	    
	    break;
	  case 5:
	    
	    break;
	  case 6:
	    secEnergyCM = recoPoint->GetEnergy6(elemId, MT, kineticE);	    
	    
	    break;
	  case 7:
	    cosLab = recoPoint->GetCos6(elemId, MT, kineticE);
	    secEnergyLab = recoPoint->GetEnergy6(elemId, MT, kineticE);	    
	    break;
	}
	break;
    }
    //double cosT = recoPoint->GetCos4(elemId, MT, kineticE);
    //double secEnergy = recoPoint->GetEnergy5(elemId, MT, kineticE);
    std::cout<<"mass = "<< particle[elemId].mass << std::endl;
    std::cout <<" MT  "<< recoPoint->MtValues[elemId][isel] <<" E= "<< kineticE <<" cos "<< cosCM 
        <<" ECM= "<< secEnergyCM<<" cosL "<< cosLab <<" Elab "<< secEnergyLab << std::endl;
    crs.shrink_to_fit();
    counter++;
  }while(counter < 1000);
  std::cout<<" end sampling "<< std::endl;
}
//------------------------------------------------------------------------------------------------------

TNudySampling::~TNudySampling(){}
