// 	This class is reconstructing probability tables for Angular distribution
// 	of the secondatries
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: March 24, 2016

#include "TList.h"
#include "Geant/TNudyEndfAng.h"
#include "Geant/TNudyEndfFile.h"
#include "Geant/TNudyEndfCont.h"
#include "Geant/TNudyEndfTab1.h"
#include "Geant/TNudyEndfTab2.h"
#include "Geant/TNudyEndfList.h"
#include "Geant/TNudyCore.h"
#include "Math/SpecFuncMathMore.h"

using namespace Nudy;
using namespace NudyPhysics;

#ifdef USE_ROOT
ClassImp(TNudyEndfAng)
#endif

TNudyEndfAng::TNudyEndfAng() : TNudyEndfSigma()
{
}
//______________________________________________________________________________
TNudyEndfAng::TNudyEndfAng(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    TIter recIter(sec->GetRecords());
    TNudyEndfCont *header = (TNudyEndfCont *)recIter.Next();
    fNLTT   = 0; 
    fLTT    = sec->GetL2();
    mtf[0]  = sec->GetMAT();
    mtf[1]  = sec->GetMT();
    mtf[2]  = sec->GetMF();
    int LI  = header->GetL1();
    switch (LI) {
      case 0:
      {
        switch (fLTT) {
          case 1: // Legendre polynomial coefficients
          {
            sec->SetCont(sec->GetC1(),sec->GetC2(),0,2,0,0);
            LegendreToCosine(recIter,sec);
          }break;
          case 2: // Tabulated probability tables
          {
            CosineToLinearCosine(recIter);
          }break;
          case 3: // Low energy given by legendre polynomial and high energy by tabulated probability tables
          {
            sec->SetCont(sec->GetC1(),sec->GetC2(),0,2,0,0);
            LegendreToCosinePlusCosine(recIter,sec);
          }break;
        }
      }break;
      case 1: // isotropic distribution
      {
        sec->SetCont(sec->GetC1(),sec->GetC2(),0,2,0,0);
        header->SetCont(header->GetC1(),header->GetC2(),0,header->GetL2(),0,0);
        // creating Tab2 parameters
        TNudyEndfTab2 *secTab2 = new TNudyEndfTab2();
        int NE = 2;
        CreateTab2(secTab2,NE);
        sec->AddAfter(header,secTab2);
        fEin.push_back(1E-5);
        fEin.push_back(1.5E8);
        for (int j = 0; j < 2; j++) {
          fCosFile4.push_back(1);
          fCosPdfFile4.push_back(0.5);
          fCosFile4.push_back(-1);
          fCosPdfFile4.push_back(0.5);
          // creating Tab1 parameters
          TNudyEndfTab1 *secTab1 = new TNudyEndfTab1();
          CreateTab1(secTab1, fCosFile4, fCosPdfFile4, mtf, 2, fEin[j]);
          sec->Add(secTab1);
          fCosFile4.clear();
          fCosPdfFile4.clear();
        }
        fEin.clear();
      }break;
    }
  } // end while loop
} // end class

TNudyEndfAng::~TNudyEndfAng()
{
  fNbt1.shrink_to_fit();
  fInt1.shrink_to_fit();
  fCosFile4.shrink_to_fit();
  fCosPdfFile4.shrink_to_fit();
  fEin.shrink_to_fit();
  fLegendCoef.shrink_to_fit();
}
//______________________________________________________________________________
void TNudyEndfAng::LegendreToCosine(TIter &recIter, TNudyEndfSec *sec)
{
  TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
  tab2->SetCont(0,0,0,0,1,tab2->GetNZ());
  tab2->SetNBT(tab2->GetNZ(),0);
  tab2->SetINT(2,0);
  fNLTT = tab2->GetNZ();
  for (int i = 0, e2N2 = tab2->GetNZ(); i != e2N2; ++i) {
    TNudyEndfList *list = (TNudyEndfList *)recIter.Next();
    fEin.push_back(list->GetC2());
    for (int j = 0, eNPL = list->GetNPL(); j != eNPL; ++j) {
      fLegendCoef.push_back(list->GetLIST(j));
    }
    int k1     = 0;
    double fme = 0.0;
    do {
      fme      = 1.0;
      double x = -1. + k1 * 0.02;
      for (int j = 0, LegendCoefSize = fLegendCoef.size(); j != LegendCoefSize; ++j) {
        double leg = ROOT::Math::legendre(j + 1, x);
        fme += 0.5 * (2. * (j + 1) + 1.) * fLegendCoef[j] * leg;
      }
      if (fme > 0.0) {
        fCosFile4.push_back(x);
        fCosPdfFile4.push_back(fme);
      }
      k1++;
    } while (k1 < 101);
    for (int l = 0; l < 100; l++) {
      RecursionLinearLeg1D(fCosFile4[l], fCosFile4[l + 1], fCosPdfFile4[l], fCosPdfFile4[l + 1]);
    }
    TNudyCore::Instance()->Sort(fCosFile4, fCosPdfFile4);
    // creating Tab1 parameters and inserting after tab1 as required by LF = 1
    TNudyEndfTab1 *secTab1 = new TNudyEndfTab1();
    CreateTab1(secTab1, fCosFile4, fCosPdfFile4, mtf, 2, fEin[i]);
    sec->AddBefore(list,secTab1);
    sec->RemoveObj(list);
    fCosFile4.clear();
    fCosPdfFile4.clear();
    fLegendCoef.clear();
  }
  fEin.clear();
}
//______________________________________________________________________________
void TNudyEndfAng::CosineToLinearCosine(TIter &recIter)
{
  TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
  for (int i = 0, e2N2 = tab2->GetNZ(); i != e2N2; ++i) {
    TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
    fEin.push_back(tab1->GetC2());
    ProcessTab1(tab1, fNR, fNP, fNbt1, fInt1, fCosFile4, fCosPdfFile4);
    int law = 2;
    for (int l = 0; l < fNP - 1; l++) {
      for (int i = 0; i < fNR; i++) {
        if (fNbt1[i] > l) {
          law = fInt1[i];
          break;
        }
      }
      if (law == 1) continue;
      RecursionLinearProb(fCosFile4[l], fCosFile4[l + 1], fCosPdfFile4[l], fCosPdfFile4[l + 1]);
    }
    TNudyCore::Instance()->Sort(fCosFile4, fCosPdfFile4);
    ModifyTab1(tab1, fCosFile4, fCosPdfFile4, fEin[i]);
    fCosFile4.clear();
    fCosPdfFile4.clear();
    fNbt1.clear();
    fInt1.clear();
  }
  fEin.clear();
}
//______________________________________________________________________________
void TNudyEndfAng::LegendreToCosinePlusCosine(TIter &recIter, TNudyEndfSec *sec)
{
  TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
  tab2->SetCont(0,0,0,0,1,tab2->GetNZ());
  tab2->SetNBT(tab2->GetNZ(),0);
  tab2->SetINT(2,0);
  fNLTT = tab2->GetNZ();
  for (int i = 0, e2N2 = tab2->GetNZ(); i != e2N2; ++i) {
    TNudyEndfList *list = (TNudyEndfList *)recIter.Next();
    fEin.push_back(list->GetC2());
    for (int j = 0, eNPL = list->GetNPL(); j != eNPL; ++j) {
      fLegendCoef.push_back(list->GetLIST(j));
    }
    int k1     = 0;
    double fme = 0.0;
    do {
      fme      = 1.0;
      double x = -1. + k1 * 0.02;
      for (int j = 0, LegendCoefSize = fLegendCoef.size(); j != LegendCoefSize; ++j) {
        double leg = ROOT::Math::legendre(j + 1, x);
        fme += 0.5 * (2. * (j + 1) + 1.) * fLegendCoef[j] * leg;
      }
      if (fme > 0.0) {
        fCosFile4.push_back(x);
        fCosPdfFile4.push_back(fme);
      }
      k1++;
    } while (k1 < 101);
    for (int l = 0; l < 100; l++) {
      RecursionLinearLeg1D(fCosFile4[l], fCosFile4[l + 1], fCosPdfFile4[l], fCosPdfFile4[l + 1]);
    }
    TNudyCore::Instance()->Sort(fCosFile4, fCosPdfFile4);
    // creating Tab1 parameters and inserting after tab1 as required by LF = 1
    TNudyEndfTab1 *secTab1 = new TNudyEndfTab1();
    CreateTab1(secTab1, fCosFile4, fCosPdfFile4, mtf, 2, fEin[i]);
    sec->AddBefore(list,secTab1);
    sec->RemoveObj(list);
    fCosFile4.clear();
    fCosPdfFile4.clear();
    fLegendCoef.clear();
  }
  TNudyEndfTab2 *tab22 = (TNudyEndfTab2 *)recIter.Next();
  int e2N2 = tab22->GetNZ();
  sec->RemoveObj(tab22);
  for (int i = 0; i != e2N2; ++i) {
    TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
    fEin.push_back(tab1->GetC2());
    ProcessTab1(tab1, fNR, fNP, fNbt1, fInt1, fCosFile4, fCosPdfFile4);
    int law = 2;
    for (int l = 0; l < fNP - 1; l++) {
      for (int i = 0; i < fNR; i++) {
        if (fNbt1[i] > l) {
          law = fInt1[i];
          break;
        }
      }
      if (law == 1) continue;
      RecursionLinearProb(fCosFile4[l], fCosFile4[l + 1], fCosPdfFile4[l], fCosPdfFile4[l + 1]);
    }
    TNudyCore::Instance()->Sort(fCosFile4, fCosPdfFile4);
    ModifyTab1(tab1, fCosFile4, fCosPdfFile4, fEin[fNLTT + i]);
    fCosFile4.clear();
    fCosPdfFile4.clear();
    fNbt1.clear();
    fInt1.clear();
  }
  tab2->SetCont(0,0,0,0,1,fNLTT + e2N2);
  tab2->SetNBT(fNLTT + e2N2,0); 
  tab2->SetINT(2,0); 
  fEin.clear();
}
//--------------------------------------------------------------------------------------
double TNudyEndfAng::RecursionLinearLeg1D(double x1, double x2, double pdf1, double pdf2)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2) return 0;
  for (int j = 0, LegendCoefSize = fLegendCoef.size(); j != LegendCoefSize; ++j) {
    double leg = ROOT::Math::legendre(j + 1, mid);
    pdf += 0.5 * (2. * (j + 1) + 1.) * fLegendCoef[j] * leg;
  }
  double pdfmid1 = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  if (pdf > 0 && fabs((pdf - pdfmid1) / pdf) <= 5E-3) {
    return 0;
  }
  fCosFile4.push_back(mid);
  fCosPdfFile4.push_back(pdf);
  RecursionLinearLeg1D(x1, mid, pdf1, pdf);
  RecursionLinearLeg1D(mid, x2, pdf, pdf2);
  return 0;
}
//--------------------------------------------------------------------------------------
double TNudyEndfAng::RecursionLinearProb(double x1, double x2, double pdf1, double pdf2)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2) return 0;
  pdf            = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, fCosFile4, fCosPdfFile4, fNP, mid);
  double pdfmid1 = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  if (pdf > 0 && fabs((pdf - pdfmid1) / pdf) <= 5E-3) {
    return 0;
  }
  fCosFile4.push_back(mid);
  fCosPdfFile4.push_back(pdf);
  RecursionLinearProb(x1, mid, pdf1, pdf);
  RecursionLinearProb(mid, x2, pdf, pdf2);
  return 0;
}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfAng::ProcessTab1(Nudy::TNudyEndfTab1 *tab1,int &NR,int &NP,rowint &fnbt, 
                               rowint &fint,rowd &x1, rowd &x2)
{  
  NudyPhysics::TNudyEndfSigma::ProcessTab1(tab1, NR, NP, fnbt, fint, x1, x2);
}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfAng::ModifyTab1(TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, double &x)
{
  NudyPhysics::TNudyEndfSigma::ModifyTab1(secTab1, x1, x2, x);
}  
// -------------------------------------------------------------------------------------------------------
void TNudyEndfAng::CreateTab2(TNudyEndfTab2 *secTab2, int &NE)
{
  NudyPhysics::TNudyEndfSigma::CreateTab2(secTab2, NE);
}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfAng::CreateTab1(TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, int mtf[], int law, double &x)
{
  NudyPhysics::TNudyEndfSigma::CreateTab1(secTab1, x1, x2, mtf, law, x);
}
