#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <iterator>
#include <algorithm>
#include <iomanip>
#include "TList.h"
#include "TFile.h"
#include "TKey.h"
#include "TTree.h"
#include "TH1D.h"
#include "Geant/TNudyEndfFile.h"
#include "Geant/TNudyEndfTape.h"
#include "Geant/TNudyEndfTab1.h"
#include "Geant/TNudyEndfTab2.h"
#include "Geant/TNudyEndfMat.h"
#include "Geant/TNudyEndfRecoPoint.h"
using namespace NudyPhysics;
int main(int, char *argv[]) {
  std::clock_t startTime = clock();
  const char* rENDF = argv[1];
  int MT, LCT;
  std::string Reaction, str, rootstr;
  double Ein, Ang, Eout;
  double Energy[5] = {1E-5, 1E3, 1E6, 1E7, 2E7};
  str = argv[1];
  std::size_t found = str.find_last_of(".");
  rootstr = str.substr(0,found)+"ThetaE.root";
  TFile *thetae = new TFile(rootstr.c_str(),"RECREATE","Angle Energy");
  TTree *ang =new TTree("ang","a Tree with angular distribution");
  ang->Branch("MT",&MT, "MT/I");  
  ang->Branch("LCT",&LCT,"LCT/I");  
  ang->Branch("Ein",&Ein,"Ein/D");  
  ang->Branch("Ang",&Ang,"Ang/D");  
  TTree *eout =new TTree("eout","a Tree with energy distribution");
  eout->Branch("MT",&MT,"MT/I");  
  eout->Branch("Ein",&Ein,"Ein/D");  
  eout->Branch("Eout",&Eout,"Eout/D");  
  TNudyEndfRecoPoint *recopoint = new TNudyEndfRecoPoint(0,rENDF);
  std::cout<<"Total XSec at 1 MeV = "<< recopoint->GetSigmaTotal(0,1E6) << "  barns" << std::endl;
  std::cout<<"Elastic XSec at 1 MeV = "<< recopoint->GetSigmaPartial(0,2,1E6) << "  barns" << std::endl;
  std::cout<<"Capture XSec at 1 MeV = "<< recopoint->GetSigmaPartial(0,102,1E6) << "  barns" << std::endl;
  for (int j = 0; j != 5; ++j) {
    for (int i = 0; i != 10000; ++i) {
      for (int mt = 0, mtSize = recopoint->fMtValues[0].size(); mt != mtSize; ++mt) {
        MT     = recopoint->fMtValues[0][mt];
        Ein    = Energy[j];
        Ang    = recopoint->GetCos4(0,MT,Ein);
        ang->Fill();
        Eout   = recopoint->GetEnergy5(0,MT,Ein);
        eout->Fill();
        // double tempE = 15E6 - Eout;
        // Eout         = tempE;
        // eout->Fill();
        //     Reaction = "Fission";
        //     Ein    = 1;
        //     Eout   = recopoint->GetEnergy5(0,18,1E6);
        //     eout->Fill();
      }
    }
  }
  thetae->Write();
  thetae->Close();
  std::cout<<"Elastic Energy is to be calculated from the kinematics "<< std::endl;
  std::clock_t endTime = clock();
  std::clock_t clockTicksTaken = endTime - startTime;
  double timeInSeconds = clockTicksTaken / (double) CLOCKS_PER_SEC;
  std::cout<<"Time in seconds = "<< timeInSeconds <<std::endl;
}
