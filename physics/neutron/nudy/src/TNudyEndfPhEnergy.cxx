//      This class is reconstructing probability tables for Energy distribution
//      of the secondatries
//      Author: Dr. Harphool Kumawat
//      Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
//      date of creation: March 24, 2016

#include "TList.h"
#include "Geant/TNudyEndfFile.h"
#include "Geant/TNudyEndfTab1.h"
#include "Geant/TNudyEndfTab2.h"
#include "Geant/TNudyCore.h"
#include "Geant/TNudyEndfPhEnergy.h"
#include "Math/SpecFuncMathMore.h"
#include "TMath.h"
using namespace Nudy;
using namespace NudyPhysics;

#ifdef USE_ROOT
ClassImp(TNudyEndfPhEnergy)
#endif

TNudyEndfPhEnergy::TNudyEndfPhEnergy() : TNudyEndfSigma()
{
}
//------------------------------------------------------------------------------------------------------
TNudyEndfPhEnergy::TNudyEndfPhEnergy(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    TIter recIter(sec->GetRecords());
    int NC                = sec->GetN1();
    for (int k = 0; k < NC; k++) {
      TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
      fNr1       = tab1->GetN1();
      fNp1       = tab1->GetN2();
      // arbitrary tabulated function (only LF=1 is defined at this moment)
      ProcessTab1(tab1, fNr1, fNp1, fNbt1, fInt1, fE1, fP1);
      TNudyEndfTab2 *tab2  = (TNudyEndfTab2 *)recIter.Next();
      ProcessTab2(tab2, fNr2, fNp2, fNbt2, fInt2);
      for (int cr = 0; cr < fNp2; cr++) {
        TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
        fEin.push_back(tab12->GetC2());
        ProcessTab1(tab12, fNr3, fNp3, fNbt3, fInt3, fEnergyFile5, fEnergyPdfFile5);
        int law = 2;
        for (int cr1 = 0; cr1 < fNp3 - 1; cr1++) {
          for (int i = 0; i < fNr3; i++) {
            if (fNbt3[i] > cr1) {
              law = fInt3[i];
              break;
            }
          }
          if (law == 1) continue;
          RecursionLinearFile5Prob(fEnergyFile5[cr1], fEnergyFile5[cr1 + 1], fEnergyPdfFile5[cr1], 
                                    fEnergyPdfFile5[cr1 + 1]);
        }
        TNudyCore::Instance()->Sort(fEnergyFile5, fEnergyPdfFile5);
        ModifyTab1(tab12, fEnergyFile5, fEnergyPdfFile5, fEin[cr]);
        fEnergyFile5.clear();
        fEnergyPdfFile5.clear();
        fNbt3.clear();
        fInt3.clear();
      }
      fNbt1.clear();
      fInt1.clear();
      fNbt2.clear();
      fInt2.clear();
      fE1.clear();
      fP1.clear();
      fEin.clear();
    }
  }
}
//------------------------------------------------------------------------------------------------------
TNudyEndfPhEnergy::~TNudyEndfPhEnergy()
{
  fE1.shrink_to_fit();
  fP1.shrink_to_fit();
  fE2.shrink_to_fit();
  fP2.shrink_to_fit();
  fE3.shrink_to_fit();
  fP3.shrink_to_fit();
  fNbt1.shrink_to_fit();
  fInt1.shrink_to_fit();
  fNbt2.shrink_to_fit();
  fInt2.shrink_to_fit();
  fNbt3.shrink_to_fit();
  fInt3.shrink_to_fit();
  fEnergyFile5.shrink_to_fit();
  fEnergyPdfFile5.shrink_to_fit();
  fEin.shrink_to_fit();
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfPhEnergy::RecursionLinearFile5Prob(double x1, double x2, double pdf1, double pdf2)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2) return 0;
  pdf            = TNudyCore::Instance()->Interpolate(fNbt3, fInt3, fNr3, fEnergyFile5, fEnergyPdfFile5, fNp3, mid);
  double pdfmid1 = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  if (fabs((pdf - pdfmid1) / pdfmid1) <= 5E-3) {
    return 0;
  }
  fEnergyFile5.push_back(mid);
  fEnergyPdfFile5.push_back(pdf);
  RecursionLinearFile5Prob(x1, mid, pdf1, pdf);
  RecursionLinearFile5Prob(mid, x2, pdf, pdf2);
  return 0;
}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfPhEnergy::ProcessTab2(Nudy::TNudyEndfTab2 *tab2, int &NR,int &NP,rowint &fnbt, rowint &fint)
{
  NudyPhysics::TNudyEndfSigma::ProcessTab2(tab2, NR, NP, fnbt, fint);
}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfPhEnergy::ProcessTab1(Nudy::TNudyEndfTab1 *tab1,int &NR,int &NP,rowint &fnbt, 
                                  rowint &fint,rowd &x1, rowd &x2)
{  
  NudyPhysics::TNudyEndfSigma::ProcessTab1(tab1, NR, NP, fnbt, fint, x1, x2);
}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfPhEnergy::ModifyTab1(TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, double &x)
{
  NudyPhysics::TNudyEndfSigma::ModifyTab1(secTab1, x1, x2, x);
}  