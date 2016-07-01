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
#include "TRandom3.h"
#endif
#include "TCanvas.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TFile.h"
#ifdef USE_ROOT
ClassImp(TNudySampling)
#endif

TNudySampling::TNudySampling(){}
//------------------------------------------------------------------------------------------------------
TNudySampling::TNudySampling(Particle* particle, TNudyEndfRecoPoint *recoPoint){
  int elemId = 0;
  double kineticE = particle[elemId].energy; 
  int isel = 0;
  std::vector<double> crs;
  int counter =0;
  double cosCM=0, cosLab=0, secEnergyCM=0, secEnergyLab=0;
  double x[100000], y[100000];
  int ecounter = 0;
  fRnd = new TRandom3(0);
  TFile *f = new TFile("test.root","recreate");
  f->cd();
  TH2D * h  = new TH2D("h", "", 1000, 1E-13, 0.2, 1000, -1, 1);
  //std::cout <<"sampling sigmaTotal "<< recoPoint->GetSigmaTotal(0,20) << std::endl;
  //std::cout <<"sampling sigmaPartial total "<< recoPoint->GetSigmaPartial(0,0,20) << std::endl;
  //std::cout <<"sampling sigmaPartial elstic "<< recoPoint->GetSigmaPartial(0,1,20) << std::endl;
  // determining reaction type from element;
  //double density = 1;
  //double charge = 1;
  //double avg = 6.022E23;
  //double ro = avg * density / mass;
  //int enemax = recoPoint->eneUni[elemId].size(); 
  do
  {
    double sum1 = 0;
    //kineticE = fRnd->Uniform(1) * recoPoint->eneUni[elemId][enemax-1];
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
    //std::cout <<" MF "<< MF << std::endl;
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
    //std::cout<<"mass = "<< particle[elemId].mass << std::endl;
    if(MT==2){
      h->Fill(secEnergyLab/1E9 , cosLab);
      x[ecounter] = secEnergyLab/1E9 ;
      y[ecounter] = cosLab ;
      ecounter ++ ;
      std::cout <<secEnergyLab/1E9 <<"  "<< cosLab <<  std::endl;
    }
    crs.clear();
    counter++;
  }while(counter < 100000);
  h->Draw("colz");
  TGraph *gr1 = new TGraph (ecounter, x , y);
  gr1->GetYaxis()->SetTitle("Cos");
  gr1->GetXaxis()->SetTitle("Energy, GeV");
  gr1->GetXaxis()->SetTitleSize(0.05);
  gr1->GetYaxis()->SetTitleSize(0.05);
  gr1->GetYaxis()->SetTitleOffset(0.5);
  gr1->GetXaxis()->SetTitleOffset(0.85);
  gr1->GetXaxis()->SetTitleSize(0.05);
  gr1->Draw();
  gr1->Write("mygraph");
  f->Write();
  std::cout<< std::endl;
}
//------------------------------------------------------------------------------------------------------

TNudySampling::~TNudySampling(){
}
