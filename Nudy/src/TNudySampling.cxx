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
  kineticE = particle[elemId].energy; 
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
    MT = recoPoint->MtValues[elemId][isel];
    MF = recoPoint->GetMt456 ( elemId, MT ) ;
    //std::cout <<" MF "<< MF <<" MT "<< MT<<"   "<< recoPoint->GetSigmaPartial(elemId,isel,kineticE) <<"  "<< kineticE << std::endl;
    // selection of the data from file 4 5 and 6 for angle and energy calculations
    //LCT = recoPoint->GetCos4Lct( elemId, MT );
    switch (MT) {
      case 2: // elastic
	residueA = particle[elemId].mass;
	residueZ = particle[elemId].charge;
	cosCM = recoPoint->GetCos4(elemId, MT, kineticE);
	cosLab = TNudyCore::Instance()->cmToLabElasticCosT(cosCM, particle[elemId].mass);
	secEnergyLab = TNudyCore::Instance()->cmToLabElasticE(kineticE, cosCM, particle[elemId].mass);
	break;
      case 11: // 2nd
	residueA = particle[elemId].mass - 3 ;
	residueZ = particle[elemId].charge - 1 ;
	break;
      case 16: //2n
	residueA = particle[elemId].mass - 1 ;
	residueZ = particle[elemId].charge ;
	GetSecParameter(particle, recoPoint);
	break;
      case 17: //3n
	residueA = particle[elemId].mass - 2 ;
	residueZ = particle[elemId].charge ;	
	GetSecParameter(particle, recoPoint);
	break;
      case 18: //fission
      {
	double nut = recoPoint->GetNuTotal(elemId, kineticE);
	int nu = (int)nut;
	if(fRnd->Uniform(1) < nut - nu) nu = nu + 1;
        for (int nui = 0; nui < nu; nui++){
	  GetSecParameter(particle, recoPoint);
	}
	double zaf = recoPoint->GetFisYield(elemId, kineticE);
      }
	break;
      case 22: //n+alpha
	residueA = particle[elemId].mass - 4 ;
	residueZ = particle[elemId].charge - 2 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 23: //n+3a
	residueA = particle[elemId].mass - 12 ;
	residueZ = particle[elemId].charge - 6 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 24: // 2n+a
	residueA = particle[elemId].mass - 5 ;
	residueZ = particle[elemId].charge - 2 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 25: //3n+a
	residueA = particle[elemId].mass - 6 ;
	residueZ = particle[elemId].charge - 2 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 28: //n+p
	residueA = particle[elemId].mass - 1 ;
	residueZ = particle[elemId].charge - 1 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 29: //n +2a
	residueA = particle[elemId].mass - 8 ;
	residueZ = particle[elemId].charge - 4 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 30: //2n+2a
	residueA = particle[elemId].mass - 9 ;
	residueZ = particle[elemId].charge - 4 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 32: //n+d
	residueA = particle[elemId].mass - 2 ;
	residueZ = particle[elemId].charge - 1 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 33: //n+t
	residueA = particle[elemId].mass - 3 ;
	residueZ = particle[elemId].charge - 1 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 34: //n + He3
	residueA = particle[elemId].mass - 3 ;
	residueZ = particle[elemId].charge - 2 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 35: //n+d+2a
	residueA = particle[elemId].mass - 10 ;
	residueZ = particle[elemId].charge - 5 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 36: //n+t+2a
	residueA = particle[elemId].mass - 11 ;
	residueZ = particle[elemId].charge - 5 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 37: //4n
	residueA = particle[elemId].mass - 3 ;
	residueZ = particle[elemId].charge ;
	GetSecParameter(particle, recoPoint);
	break;
      case 41: //2n+p
	residueA = particle[elemId].mass - 2 ;
	residueZ = particle[elemId].charge - 1 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 42: //3n+p
	residueA = particle[elemId].mass - 3 ;
	residueZ = particle[elemId].charge - 1 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 44: //n+2p
	residueA = particle[elemId].mass - 2 ;
	residueZ = particle[elemId].charge - 1 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 45: //n+p+a
	residueA = particle[elemId].mass - 5 ;
	residueZ = particle[elemId].charge - 3 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 50:
      case 51:
      case 52:
      case 53:
      case 54:
      case 55:
      case 56:
      case 57:
      case 58:
      case 59:
      case 60:
      case 61:
      case 62:
      case 63:
      case 64:
      case 65:
      case 66:
      case 67:
      case 68:
      case 69:
      case 70:
      case 71:
      case 72:
      case 73:
      case 74:
      case 75:
      case 76:
      case 77:
      case 78:
      case 79:
      case 80:
      case 81:
      case 82:
      case 83:
      case 84:
      case 85:
      case 86:
      case 87:
      case 88:
      case 89:
      case 90:
      case 91:
	residueA = particle[elemId].mass ;
	residueZ = particle[elemId].charge ;
	cosCM = recoPoint->GetCos4(elemId, MT, kineticE);
	secEnergyCM = recoPoint->GetEnergy5(elemId, MT, kineticE);
	secEnergyLab = TNudyCore::Instance()->cmToLabInelasticE(secEnergyCM, kineticE, cosCM, particle[elemId].mass);
	cosLab = TNudyCore::Instance()->cmToLabInelasticCosT(secEnergyLab, secEnergyCM, kineticE, cosCM, particle[elemId].mass);	  
	break;
      case 102: //capture
	residueA = particle[elemId].mass + 1 ;
	residueZ = particle[elemId].charge ;
	break;
      case 103: //p
	residueA = particle[elemId].mass ;
	residueZ = particle[elemId].charge - 1 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 104: //d
	residueA = particle[elemId].mass - 1 ;
	residueZ = particle[elemId].charge - 1 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 105: //t
	residueA = particle[elemId].mass - 2 ;
	residueZ = particle[elemId].charge - 1 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 106: //He3
	residueA = particle[elemId].mass - 2 ;
	residueZ = particle[elemId].charge - 2 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 107: //alpha
	residueA = particle[elemId].mass - 3 ;
	residueZ = particle[elemId].charge - 2 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 108: //2a
	residueA = particle[elemId].mass - 7 ;
	residueZ = particle[elemId].charge - 4 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 109: //3a
	residueA = particle[elemId].mass - 11 ;
	residueZ = particle[elemId].charge - 6 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 111: //2p
	residueA = particle[elemId].mass - 1 ;
	residueZ = particle[elemId].charge - 2 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 112: //p+a
	residueA = particle[elemId].mass - 4 ;
	residueZ = particle[elemId].charge - 3 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 113: //t+2a
	residueA = particle[elemId].mass - 10 ;
	residueZ = particle[elemId].charge - 5 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 114: //d+2a
	residueA = particle[elemId].mass - 9 ;
	residueZ = particle[elemId].charge - 5 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 115: //p+d
	residueA = particle[elemId].mass - 2 ;
	residueZ = particle[elemId].charge - 2 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 116: //p+t
	residueA = particle[elemId].mass - 3 ;
	residueZ = particle[elemId].charge - 2 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 117: //d+a
	residueA = particle[elemId].mass - 5 ;
	residueZ = particle[elemId].charge - 3 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 152: //5n
	residueA = particle[elemId].mass - 4 ;
	residueZ = particle[elemId].charge ;
	GetSecParameter(particle, recoPoint);
	break;
      case 153: //6n
	residueA = particle[elemId].mass - 5 ;
	residueZ = particle[elemId].charge ;
	GetSecParameter(particle, recoPoint);
	break;
      case 154: //2n+t
	residueA = particle[elemId].mass - 4 ;
	residueZ = particle[elemId].charge - 1 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 155: //t+a
	residueA = particle[elemId].mass - 6 ;
	residueZ = particle[elemId].charge - 3 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 156: //4n+p
	residueA = particle[elemId].mass - 4 ;
	residueZ = particle[elemId].charge - 1 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 157: //3n+d
	residueA = particle[elemId].mass - 4 ;
	residueZ = particle[elemId].charge - 1 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 158: //n+d+a
	residueA = particle[elemId].mass - 6 ;
	residueZ = particle[elemId].charge - 3 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 159: //2n+p+a
	residueA = particle[elemId].mass - 6 ;
	residueZ = particle[elemId].charge - 3 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 160: //7n
	residueA = particle[elemId].mass - 6 ;
	residueZ = particle[elemId].charge ;
	GetSecParameter(particle, recoPoint);
	break;
      case 161: //8n
	residueA = particle[elemId].mass - 7 ;
	residueZ = particle[elemId].charge ;
	GetSecParameter(particle, recoPoint);
	break;
      case 162: //5np
	residueA = particle[elemId].mass - 5 ;
	residueZ = particle[elemId].charge - 1 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 163: //6np
	residueA = particle[elemId].mass - 6 ;
	residueZ = particle[elemId].charge - 1 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 164: //7np
	residueA = particle[elemId].mass - 7 ;
	residueZ = particle[elemId].charge - 1 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 165: //4n+a
	residueA = particle[elemId].mass - 7 ;
	residueZ = particle[elemId].charge - 2 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 166: //5na
	residueA = particle[elemId].mass - 8 ;
	residueZ = particle[elemId].charge - 2 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 167: //6na
	residueA = particle[elemId].mass - 9 ;
	residueZ = particle[elemId].charge - 2 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 168: //7na
	residueA = particle[elemId].mass - 10 ;
	residueZ = particle[elemId].charge - 2 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 169: //4nd
	residueA = particle[elemId].mass - 5 ;
	residueZ = particle[elemId].charge - 1 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 170: //5nd
	residueA = particle[elemId].mass - 6 ;
	residueZ = particle[elemId].charge - 1 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 171: //6nd
	residueA = particle[elemId].mass - 7 ;
	residueZ = particle[elemId].charge - 1 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 172: //3nt
	residueA = particle[elemId].mass - 5 ;
	residueZ = particle[elemId].charge - 1 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 173: //4nt
	residueA = particle[elemId].mass - 6 ;
	residueZ = particle[elemId].charge - 1 ;
	GetSecParameter(particle, recoPoint);
	break; 
      case 174: //5nt
	residueA = particle[elemId].mass - 7 ;
	residueZ = particle[elemId].charge - 1 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 175: //6nt
	residueA = particle[elemId].mass - 8 ;
	residueZ = particle[elemId].charge - 1 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 176: //2n+He3
	residueA = particle[elemId].mass - 4 ;
	residueZ = particle[elemId].charge - 2 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 177: //3n + He3
	residueA = particle[elemId].mass - 5 ;
	residueZ = particle[elemId].charge - 2 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 178: //4n +He3
	residueA = particle[elemId].mass - 6 ;
	residueZ = particle[elemId].charge - 2 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 179: //3n2p
	residueA = particle[elemId].mass - 4 ;
	residueZ = particle[elemId].charge - 2 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 180: //3n2a
	residueA = particle[elemId].mass - 10 ;
	residueZ = particle[elemId].charge - 4 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 181: //3npa
	residueA = particle[elemId].mass - 7 ;
	residueZ = particle[elemId].charge - 3 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 182: //dt
	residueA = particle[elemId].mass - 4 ;
	residueZ = particle[elemId].charge - 2 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 183: //npd
	residueA = particle[elemId].mass - 3 ;
	residueZ = particle[elemId].charge - 2 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 184: //npt
	residueA = particle[elemId].mass - 4 ;
	residueZ = particle[elemId].charge - 2 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 185: //ndt
	residueA = particle[elemId].mass - 5 ;
	residueZ = particle[elemId].charge - 2 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 186: //npHe3
	residueA = particle[elemId].mass - 4 ;
	residueZ = particle[elemId].charge - 3 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 187: //ndHe3
	residueA = particle[elemId].mass - 5 ;
	residueZ = particle[elemId].charge - 3 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 188: //ntHe3
	residueA = particle[elemId].mass - 6 ;
	residueZ = particle[elemId].charge - 3 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 189: //nta
	residueA = particle[elemId].mass - 7 ;
	residueZ = particle[elemId].charge - 3 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 190: //2n2p
	residueA = particle[elemId].mass - 3 ;
	residueZ = particle[elemId].charge - 2 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 191: //pHe3
	residueA = particle[elemId].mass - 3 ;
	residueZ = particle[elemId].charge - 3 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 192: //dHe3
	residueA = particle[elemId].mass - 4 ;
	residueZ = particle[elemId].charge - 3 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 193: //aHe3
	residueA = particle[elemId].mass - 6 ;
	residueZ = particle[elemId].charge - 4 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 194: //4n2p
	residueA = particle[elemId].mass - 5 ;
	residueZ = particle[elemId].charge - 2 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 195: //4n2a
	residueA = particle[elemId].mass - 11 ;
	residueZ = particle[elemId].charge - 4 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 196: //4npa
	residueA = particle[elemId].mass - 8 ;
	residueZ = particle[elemId].charge - 3 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 197: //3p
	residueA = particle[elemId].mass - 2 ;
	residueZ = particle[elemId].charge - 3 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 198: //n3p
	residueA = particle[elemId].mass - 3 ;
	residueZ = particle[elemId].charge - 3 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 199: //3n2pa
	residueA = particle[elemId].mass - 8 ;
	residueZ = particle[elemId].charge - 4 ;
	GetSecParameter(particle, recoPoint);
	break;
      case 200: //5n2p
	residueA = particle[elemId].mass - 6 ;
	residueZ = particle[elemId].charge - 2 ;
	GetSecParameter(particle, recoPoint);
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
void TNudySampling::GetSecParameter(Particle* particle, TNudyEndfRecoPoint *recoPoint){
  //std::cout<<"MF "<< MF << " MT "<< MT << " ID "<< elemId << std::endl;
  switch (MF) {
    case 4:
      cosCM = recoPoint->GetCos4(elemId, MT, kineticE);
      secEnergyCM = recoPoint->GetEnergy5(elemId, MT, kineticE);
      secEnergyLab = TNudyCore::Instance()->cmToLabInelasticE(secEnergyCM, kineticE, cosCM, particle[elemId].mass);
      cosLab = TNudyCore::Instance()->cmToLabInelasticCosT(secEnergyLab, secEnergyCM, kineticE, cosCM, particle[elemId].mass);	  
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
}
//__________________________________________________________________________________________________________________
TNudySampling::~TNudySampling(){
}
