//      This class is reconstructing photon multiplicities and transition arrays
//      Author: Dr. Harphool Kumawat
//      Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
//      date of creation: March 24, 2016

#include "TList.h"
#include "Geant/TNudyEndfFile.h"
#include "Geant/TNudyEndfList.h"
#include "Geant/TNudyEndfTab1.h"
#include "Geant/TNudyCore.h"
#include "Geant/TNudyEndfPhYield.h"

using namespace Nudy;
using namespace NudyPhysics;

#ifdef USE_ROOT
ClassImp(TNudyEndfPhYield)
#endif

TNudyEndfPhYield::TNudyEndfPhYield() : TNudyEndfSigma()
{
}

//______________________________________________________________________________
TNudyEndfPhYield::TNudyEndfPhYield(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    TIter recIter(sec->GetRecords());
    int NK  = sec->GetN1();
    int LO  = sec->GetL1();
    int law = 2;
    switch (LO) {
      case 1: // multiplicities
      {
        switch (NK) {
          case 1:
          {
            TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
            fEG = tab1->GetC1();
            fES = tab1->GetC2(); 
            fLP = tab1->GetL1(); 
            fLF = tab1->GetL2(); 
            ProcessTab1(tab1, fNR, fNP, fNBT, fINT, fEIN, fMultiPh);
            for (int i = 0; i != fNP -1; ++i) {
              for (int j = 0; j < fNR; j++) {
                if (fNBT[j] > i) {
                  law = fINT[j];
                  break;
                }
              }
              if (law == 1) continue;
              RecursionLinear(fEIN[i],fEIN[i + 1],fMultiPh[i],fMultiPh[i + 1]);
            }
            TNudyCore::Instance()->Sort(fEIN, fMultiPh);
            ModifyTab1(tab1, fEIN, fMultiPh, fES);
            fNP = fEIN.size();
            tab1->SetCont(fEG, fES, fLP, fLF, 1, fNP);
            tab1->SetNBT(fNP,0);
            tab1->SetINT(2,0);
            for (int crs = 0; crs < fNP; crs++) {
              tab1->SetX(fEIN[crs],crs);
              tab1->SetY(fMultiPh[crs],crs);
            }
            fEIN.clear();
            fMultiPh.clear();
          }break;
          default:
          {
            TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
            ProcessTab1(tab1, fNR, fNP, fNBT, fINT, fEIN, fMultiPh);
            fES = 0;
            for (int i = 0; i != fNP -1; ++i){
              for (int j = 0; j < fNR; j++) {
                if (fNBT[j] > i) {
                  law = fINT[j];
                  break;
                }
              }
              if (law == 1) continue;
              RecursionLinear(fEIN[i],fEIN[i + 1],fMultiPh[i],fMultiPh[i + 1]);
            }
            TNudyCore::Instance()->Sort(fEIN, fMultiPh);
            ModifyTab1(tab1, fEIN, fMultiPh, fES);
            fEIN.clear();
            fMultiPh.clear();
            for (int nki = 0; nki < NK; nki++) {
              TNudyEndfTab1 *tab11 = (TNudyEndfTab1 *)recIter.Next();
              fEG = tab11->GetC1();
              fES = tab11->GetC2();
              fLP = tab11->GetL1();
              fLF = tab11->GetL2();
              fNR = tab11->GetN1();
              fNP = tab11->GetN2();
              ProcessTab1(tab11, fNR, fNP, fNBT, fINT, fEIN, fMultiPh);
              for (int i = 0; i != fNP -1; ++i){
                for (int j = 0; j < fNR; j++) {
                  if (fNBT[j] > i) {
                    law = fINT[j];
                    break;
                  }
                }
                if (law == 1) continue;
                RecursionLinear(fEIN[i],fEIN[i + 1],fMultiPh[i],fMultiPh[i + 1]);
              }
              TNudyCore::Instance()->Sort(fEIN, fMultiPh);
              ModifyTab1(tab11, fEIN, fMultiPh, fES);
              int NP = fEIN.size();
              tab11->SetCont(fEG, fES, fLP, fLF, 1, NP);
              tab11->SetNBT(fNP,0);
              tab11->SetINT(2,0);
              for (int crs = 0; crs < fNP; crs++) {
                tab11->SetX(fEIN[crs],crs);
                tab11->SetY(fMultiPh[crs],crs);
              }
              fEIN.clear();
              fMultiPh.clear();
            }
          }break;
        }
      }break;
      case 2: // Transition Probability Arrays
      {
        int LG = sec->GetL2();
        TNudyEndfList *list1 = (TNudyEndfList *)recIter.Next();
        int NT = list1->GetN2();
        switch (LG) { 
          case 1: // all transitions are gamma- transitions only
          {
            for (int j = 0; j < NT; j++) {
              fENS.push_back(list1->GetLIST(j * 2));
              fTP.push_back(list1->GetLIST(j * 2 + 1));
            } 
            fENS.clear();
            fTP.clear();
          }break;
          case 2: // complex case where internal transition also occur
          {
            for (int j = 0; j < NT; j++) {
              fENS.push_back(list1->GetLIST(j * 3));
              fTP.push_back(list1->GetLIST(j * 3 + 1));
              fGP.push_back(list1->GetLIST(j * 3 + 2));
            }
            fENS.clear();
            fTP.clear();
            fGP.clear();
          }break;
        }
      }break;
    }
  }
}
//------------------------------------------------------------------------------------------------------
TNudyEndfPhYield::~TNudyEndfPhYield() {}
//------------------------------------------------------------------------------------------------------
double TNudyEndfPhYield::RecursionLinear(double x1, double x2, double pdf1, double pdf2)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2) return 0;
  pdf            = TNudyCore::Instance()->Interpolate(fNBT, fINT, fNR, fEIN, fMultiPh, fNP, mid);
  double pdfmid1 = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  if (fabs((pdf - pdfmid1) / pdfmid1) <= 5E-3) {
    return 0;
  }
  fEIN.push_back(mid);
  fMultiPh.push_back(pdf);
  RecursionLinear(x1, mid, pdf1, pdf);
  RecursionLinear(mid, x2, pdf, pdf2);
  return 0;
}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfPhYield::ProcessTab1(Nudy::TNudyEndfTab1 *tab1,int &NR,int &NP,rowint &fnbt, 
                                  rowint &fint,rowd &x1, rowd &x2)
{  
  NudyPhysics::TNudyEndfSigma::ProcessTab1(tab1, NR, NP, fnbt, fint, x1, x2);
}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfPhYield::ModifyTab1(TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, double &x)
{
  NudyPhysics::TNudyEndfSigma::ModifyTab1(secTab1, x1, x2, x);
}  
