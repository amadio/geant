// 	This class is reconstructing photon production cross-section
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: March 24, 2016

#include "TList.h"
#include "Geant/TNudyEndfFile.h"
#include "Geant/TNudyEndfList.h"
#include "Geant/TNudyEndfTab1.h"
#include "Geant/TNudyCore.h"
#include "Geant/TNudyEndfPhProd.h"
using namespace Nudy;
using namespace NudyPhysics;

#ifdef USE_ROOT
ClassImp(TNudyEndfPhProd)
#endif
TNudyEndfPhProd::TNudyEndfPhProd() : TNudyEndfSigma()
{
}
//______________________________________________________________________________
TNudyEndfPhProd::TNudyEndfPhProd(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    TIter recIter(sec->GetRecords());
    int NK     = sec->GetN1();
    switch (NK) {
      case 1:
      {
        TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
        fEG = tab1->GetC1();
        fES = tab1->GetC2(); 
        fLP = tab1->GetL1(); 
        fLF = tab1->GetL2(); 
        ProcessTab1(tab1, fNR, fNP, fNBT, fINT, fEIN, fSigPh);
        int law = 2;
        for (int i = 0; i != fNP -1; ++i){
          for (int j = 0; j < fNR; j++) {
            if (fNBT[j] > i) {
              law = fINT[j];
              break;
            }
          }
          if (law == 1) continue;
          RecursionLinear(fEIN[i], fEIN[i + 1], fSigPh[i], fSigPh[i + 1]);
        }
        TNudyCore::Instance()->Sort(fEIN, fSigPh);
        ModifyTab1(tab1, fEIN, fSigPh, fES);
        fNP = fEIN.size();
        tab1->SetCont(fEG, fES, fLP, fLF, 1, fNP);
        tab1->SetNBT(fNP,0);
        tab1->SetINT(2,0);
        for (int crs = 0; crs < fNP; crs++) {
          tab1->SetX(fEIN[crs],crs);
          tab1->SetY(fSigPh[crs],crs);
        }
        fEIN.clear();
        fSigPh.clear();
      }break;
      default:
      {
        TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
        ProcessTab1(tab1, fNR, fNP, fNBT, fINT, fEIN, fSigPh);
        fES = 0;
        for (int i = 0; i != fNP -1; ++i){
          RecursionLinear(fEIN[i], fEIN[i + 1], fSigPh[i], fSigPh[i + 1]);
        }
        TNudyCore::Instance()->Sort(fEIN, fSigPh);
        ModifyTab1(tab1, fEIN, fSigPh, fES);
        fEIN.clear();
        fSigPh.clear();
        for (int nki = 0; nki < NK; nki++) {
          TNudyEndfTab1 *tab11 = (TNudyEndfTab1 *)recIter.Next();
          fEG = tab11->GetC1();
          fES = tab11->GetC2();
          fLP = tab11->GetL1();
          fLF = tab11->GetL2();
          fNR = tab11->GetN1();
          fNP = tab11->GetN2();
          ProcessTab1(tab1, fNR, fNP, fNBT, fINT, fEIN, fSigPh);
          int law = 2;
          for (int i = 0; i != fNP -1; ++i){
            for (int j = 0; j < fNR; j++) {
              if (fNBT[j] > i) {
                law = fINT[j];
                break;
              }
            }
            if (law == 1) continue;
            RecursionLinear(fEIN[i], fEIN[i + 1], fSigPh[i], fSigPh[i + 1]);
          }
          TNudyCore::Instance()->Sort(fEIN, fSigPh);
          ModifyTab1(tab1, fEIN, fSigPh, fES);
          fNP = fEIN.size();
          tab1->SetCont(fEG, fES, fLP, fLF, 1, fNP);
          tab1->SetNBT(fNP,0);
          tab1->SetINT(2,0);
          for (int crs = 0; crs < fNP; crs++) {
            tab1->SetX(fEIN[crs],crs);
            tab1->SetY(fSigPh[crs],crs);
          }
          fEIN.clear();
          fSigPh.clear();
        }
      }break;
    }
  }
}
//------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------
double TNudyEndfPhProd::RecursionLinear(double x1, double x2, double pdf1, double pdf2)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2) return 0;
  pdf            = TNudyCore::Instance()->Interpolate(fNBT, fINT, fNR, fEIN, fSigPh, fNP, mid);
  double pdfmid1 = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  if (fabs((pdf - pdfmid1) / pdfmid1) <= 5E-3) {
    return 0;
  }
  fEIN.push_back(mid);
  fSigPh.push_back(pdf);
  RecursionLinear(x1, mid, pdf1, pdf);
  RecursionLinear(mid, x2, pdf, pdf2);
  return 0;
}
// -------------------------------------------------------------------------------------------------------
TNudyEndfPhProd::~TNudyEndfPhProd() {}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfPhProd::ProcessTab1(Nudy::TNudyEndfTab1 *tab1,int &NR,int &NP,rowint &fnbt, 
                               rowint &fint,rowd &x1, rowd &x2)
{  
  NudyPhysics::TNudyEndfSigma::ProcessTab1(tab1, NR, NP, fnbt, fint, x1, x2);
}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfPhProd::ModifyTab1(TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, double &x)
{
  NudyPhysics::TNudyEndfSigma::ModifyTab1(secTab1, x1, x2, x);
}  
