#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
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

using namespace Nudy;
using namespace NudyPhysics;
int main(int, char *argv[]) {
  std::clock_t startTime = clock();
  const char* rENDF = argv[1];
  std::string str, rootstr;
  int MT, ZAP, LCT;
  double Energy, Xsect, Ein, Ang, AngPdf, Eout, EoutPdf;
  std::vector<double> EnergyMt,XsecFiss, XsecCap;
  std::vector<double> x1,x2;
  str = argv[1];
  std::size_t found = str.find_last_of(".");
  rootstr = str.substr(0,found)+"Xsec.root";
  TFile *Xsec = new TFile(rootstr.c_str(),"RECREATE","Xsection");
  TTree *tree = new TTree("XsecTree","A ENDF ROOT tree");
  tree->Branch("MT",&MT,"MT/I");  
  tree->Branch("E",&Energy,"Energy/D");  
  tree->Branch("Xsec",&Xsect,"Xsect/D");  
  TTree *eta =new TTree("eta","a Tree with fast fission data");
  TH1D *NuTotal = new TH1D("NuTotal", "", 150, -5, 8);
  TH1D *EnergyNu = new TH1D("EnergyNu", "", 150, -5, 8);
  TH1D *FissionNuTotal = new TH1D("FissionNuTotal", "", 150, -5, 8);
  TTree *ang =new TTree("ang","a Tree with angular distribution");
  ang->Branch("MT",&MT,"MT/I");  
  ang->Branch("LCT",&LCT,"LCT/I");  
  ang->Branch("Ein",&Ein,"Ein/D");  
  ang->Branch("Ang",&Ang,"Ang/D");  
  ang->Branch("AngPdf",&AngPdf,"AngPdf/D");  
  TTree *eout =new TTree("eout","a Tree with energy distribution");
  eout->Branch("MT",&MT,"MT/I");  
  eout->Branch("LCT",&LCT,"LCT/I");  
  eout->Branch("Ein",&Ein,"Ein/D");  
  eout->Branch("Eout",&Eout,"Eout/D");  
  eout->Branch("EoutPdf",&EoutPdf,"EoutPdf/D");  
  TTree *angeout =new TTree("angeout","a Tree with Angle-energy distribution");
  angeout->Branch("MT",&MT,"MT/I");  
  angeout->Branch("LCT",&LCT,"LCT/I");  
  angeout->Branch("Ein",&Ein,"Ein/D");  
  angeout->Branch("Ang",&Ang,"Ang/D");  
  angeout->Branch("AngPdf",&AngPdf,"AngPdf/D");  
  angeout->Branch("Eout",&Eout,"Eout/D");  
  angeout->Branch("EoutPdf",&EoutPdf,"EoutPdf/D");  
  
  TFile *rEND = TFile::Open(rENDF);
  if (!rEND || rEND->IsZombie()) printf("Error: TFile :: Cannot open file %s\n", rENDF);
  TKey *rkey          = (TKey *)rEND->GetListOfKeys()->First();
  TNudyEndfTape *tape = (TNudyEndfTape *)rkey->ReadObj();
  TNudyEndfMat *tMat  = nullptr;
  TList *mats         = (TList *)tape->GetMats();
  int nmats           = mats->GetEntries();
  for (int iMat = 0; iMat < nmats; iMat++) {
    tMat = (TNudyEndfMat *)mats->At(iMat);
    TIter iter(tMat->GetFiles());
    TNudyEndfFile *file = nullptr;
    while ((file = (TNudyEndfFile *)iter.Next())) {
      switch (file->GetMF()) {
        case 1:
        {
          TIter secIter(file->GetSections());
          TNudyEndfSec *sec = nullptr;
          while ((sec = (TNudyEndfSec *)secIter.Next())) {
            MT        = sec->GetMT();
            TIter recIter(sec->GetRecords());
            if(MT !=452 )continue;
            TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
            int NP = tab1->GetNP();
            for (int crs = 0; crs < NP; crs++) {
              x1.push_back(tab1->GetX(crs));
              x2.push_back(tab1->GetY(crs));
            }
          }
        }break;
        case 3:
        {
          TIter secIter(file->GetSections());
          TNudyEndfSec *sec = nullptr;
          while ((sec = (TNudyEndfSec *)secIter.Next())) {
            MT        = sec->GetMT();
            TIter recIter(sec->GetRecords());
            TNudyEndfCont *cont = (TNudyEndfCont *)recIter.Next();
            TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)(sec->GetRecords()->At(0));
            for (int crs  = 0, NP = cont->GetN2(); crs != NP; ++crs) {
              Energy = tab1->GetX(crs);
              Xsect  = tab1->GetY(crs);
              if (MT == 18) {
                EnergyMt.push_back(Energy);
                XsecFiss.push_back(Xsect);
              } else if (MT == 102) {
                XsecCap.push_back(Xsect);
              }
              tree->Fill();
            }
          } 
        }break;
        case 6:
        {
          TIter secIter(file->GetSections());
          TNudyEndfSec *sec;
          while ((sec = (TNudyEndfSec *)secIter.Next())) {
            int NK = sec->GetN1();
            TIter recIter(sec->GetRecords());
            for (int k = 0; k < NK; k++) {
              TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
              MT      = sec->GetMT();
              ZAP     = tab1->GetC1();
              LCT     = sec->GetL2();
              int LAW     = tab1->GetL2();
              // flag to identify if it is angular from file (4, 14) or energy distribution from file (5, 15)
              int CosOrEn = sec->GetN2(); 
              // There is no structure for these laws
              if (LAW == 3 || LAW == 4 || LAW == 0) continue;
              if (ZAP != 1) continue;
              switch (LAW) {
                case 2: // this law testedwith n-005_B_010.endf file
                { 
                  switch (CosOrEn) {
                    case 0: // from file 6
                    case 4: // from file 4 
                    {
                      TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
                      int Np2   = tab2->GetNZ();
                      for (int lis = 0; lis < Np2; lis++) {
                        TNudyEndfTab1 *tab11 = (TNudyEndfTab1 *)recIter.Next();
                        Ein = tab11->GetC2();
                        for (int crs = 0; crs < tab11->GetNP(); crs++) {
                          Ang    = tab11->GetX(crs);
                          AngPdf = tab11->GetY(crs);
                          ang->Fill();  
                        }
                      }
                    }break;
                    case 5:
                    {
                      TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
                      for (int cr = 0; cr < tab2->GetNZ(); cr++) {
                        TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
                        Ein = tab12->GetC2();
                        for (int crs = 0; crs < tab12->GetNP(); crs++) {
                          Eout    = tab12->GetX(crs);
                          EoutPdf = tab12->GetY(crs);
                          eout->Fill();  
                        }
                      }
                    }break;
                  }
                }break;
                case 7: // laboratory angle and energy law
                {
                  TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
                  for (int cr1 = 0; cr1 < tab2->GetNZ(); cr1++) {
                    TNudyEndfTab2 *tab3 = (TNudyEndfTab2 *)recIter.Next();
                    Ein = tab3->GetC2();
                    int NMU = tab3->GetN2();
                    for (int cr2 = 0; cr2 < NMU; cr2++) {
                      TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
                      Ang = tab12->GetC2();
                      double fsum = 0.0;
                      for (int crs = 0; crs < tab12->GetNP(); crs++) {
                        Eout    = tab12->GetX(crs);
                        EoutPdf = tab12->GetY(crs);
                        angeout->Fill();
                        if (crs > 1) {
                          fsum += 0.5 * (tab12->GetY(crs) + tab12->GetY(crs-1))*
                                  (tab12->GetX(crs) - tab12->GetX(crs-1));
                        }
                      }
                      AngPdf = fsum;
                      angeout->Fill();
                    }
                  }
                } break;
              }
            }
          }
        }break;
      }
    }
  }  
  double sigfis, sigcap;
  double energ, etavalue;
  eta->Branch("energy",&energ,"energ/D");
  eta->Branch("eta",&etavalue,"etavalue/D");
  if (EnergyMt.size() > 1) {
    for (int j1 = 0, EnergyMtSize = EnergyMt.size(); j1 != EnergyMtSize ; ++j1){
      sigfis = XsecFiss[j1];
      sigcap = XsecCap[j1];
      int min = 0;
      int max = x1.size() - 1;
      int mid = 0;
      if (EnergyMt[j1] <= x1[min])
        min = 0;
      else if (EnergyMt[j1] >= x1[max])
        min = max - 1;
      else {
        while (max - min > 1) {
          mid = (min + max) / 2;
          if (EnergyMt[j1] < x1[mid])
            max = mid;
          else
            min = mid;
        }
      }
      double nutotal = x2[min] +(x2[min + 1] - x2[min]) * (EnergyMt[j1] - x1[min]) /(x1[min + 1] - x1[min]);
      if (nutotal/(1+sigcap/sigfis)>0){
        NuTotal->Fill(log10(EnergyMt[j1]),nutotal/(1+sigcap/sigfis));
        EnergyNu->Fill(log10(EnergyMt[j1]));
      }
      energ    = EnergyMt[j1];
      etavalue = nutotal/(1+sigcap/sigfis);
      for(int j = 0; j < 150; j++){
        double num = NuTotal->GetBinContent(j);
        int   num2 = EnergyNu->GetBinContent(j);
        if ( num >0 && num2 > 0)FissionNuTotal->SetBinContent(j,num/num2);
      }
      eta->Fill();
    }
  }        
  rEND->Close();
  Xsec->Write();
  Xsec->Close();
  std::clock_t endTime = clock();
  std::clock_t clockTicksTaken = endTime - startTime;
  double timeInSeconds = clockTicksTaken / (double) CLOCKS_PER_SEC;
  std::cout<<"Time in seconds = "<< timeInSeconds <<std::endl;
}
