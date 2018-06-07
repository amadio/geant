//      This class is reconstructing probability tables for Angular distribution
//      of the secondatries
//      Author: Dr. Harphool Kumawat
//      Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
//      date of creation: March 24, 2016

#include "TList.h"
#include "Geant/TNudyEndfFile.h"
#include "Geant/TNudyEndfCont.h"
#include "Geant/TNudyEndfTab1.h"
#include "Geant/TNudyEndfTab2.h"
#include "Geant/TNudyEndfList.h"
#include "Geant/TNudyCore.h"
#include "Geant/TNudyEndfMat.h"
#include "Math/SpecFuncMathMore.h"
#include "Geant/TNudyEndfPhAng.h"


using namespace Nudy;
using namespace NudyPhysics;


#ifdef USE_ROOT
ClassImp(TNudyEndfPhAng)
#endif

TNudyEndfPhAng::TNudyEndfPhAng():TNudyEndfSigma()
{
}
//______________________________________________________________________________
TNudyEndfPhAng::TNudyEndfPhAng(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    TIter recIter(sec->GetRecords());
    int LTT     = sec->GetL2();
    int LI      = sec->GetL1();
    double EGk  = 0;
    fNK         = sec->GetN1();
    fNI         = sec->GetN2();
    mtf[0]      = sec->GetMAT();
    mtf[1]      = sec->GetMT();
    mtf[2]      = sec->GetMF();
    switch (LI) {
      case 1: // isotropic distribution
      { 
        sec->SetCont(sec->GetC1(),sec->GetC2(),0,LTT,fNK,0);
        TNudyEndfTab2 *secTab2 = new TNudyEndfTab2();
        CreateTab2(secTab2,fNK);
        sec->Add(secTab2);
        TNudyEndfTab1 *secTab1(new TNudyEndfTab1[fNK]);
        for (int j = 0; j < fNK; j++) {
          fCosFile4.push_back(1);
          fCosPdfFile4.push_back(0.5);
          fCosFile4.push_back(-1);
          fCosPdfFile4.push_back(0.5);
          CreateTab1(&secTab1[j], fCosFile4, fCosPdfFile4, mtf, 2, EGk);
          if (j == 0) {
            sec->AddAfter(secTab2, &secTab1[j]);
          } else {
            sec->AddAfter(&secTab1[j-1], &secTab1[j]);
          }
          fCosFile4.clear();
          fCosPdfFile4.clear();
        }
      } break;
      case 0: // isotropic discrete and anisotropic continuous distributions
      {
        if (fNI > 0) {
          sec->SetCont(sec->GetC1(),sec->GetC2(),0,LTT,fNK,0);
          TNudyEndfTab2 *secTab2 = new TNudyEndfTab2();
          CreateTab2(secTab2,fNK);
          TNudyEndfTab1 *secTab1(new TNudyEndfTab1[fNI]);
          for (int k = 0; k != fNI; ++k) { // discrete isotropic spectra
            TNudyEndfCont *Cont1 = (TNudyEndfCont *)recIter.Next();
            if (k == 0) sec->AddBefore(Cont1, secTab2);
            double energy = Cont1->GetC1();
            fCosFile4.push_back(1);
            fCosPdfFile4.push_back(0.5);
            fCosFile4.push_back(-1);
            fCosPdfFile4.push_back(0.5);
            CreateTab1(&secTab1[k], fCosFile4, fCosPdfFile4, mtf, 2, energy);
            if (k == 0) {
              sec->AddAfter(Cont1, &secTab1[k]);
            } else {
              sec->AddAfter(&secTab1[k-1], &secTab1[k]);
            }
            fCosFile4.clear();
            fCosPdfFile4.clear();
            sec->RemoveObj(Cont1);
          }
        }
        switch (LTT) {
          case 1: // Legendre polynomial coefficients
          { 
            sec->SetCont(sec->GetC1(),sec->GetC2(),0,2,fNK,0);
            // continuous anisotropic spectra
            TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
            // NE = eN2 are the number of continuous photon energy distributions
            for (int i = 0, eN2 = tab2->GetN2(); i != eN2; ++i) {
              TNudyEndfList *list = (TNudyEndfList *)recIter.Next();
              fEin.push_back(list->GetC2());
              // NL = eNPL no. of legendre coefficients
              for (int j = 0, eNPL = list->GetNPL(); j != eNPL; ++j) {
                fLegendCoef.push_back(list->GetLIST(j));
              }
              int k1     = 0;
              double fme = 0.0;
              do {
                fme      = 0.5;
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
            if (fNI > 0) sec->RemoveObj(tab2);
          } break;
          case 2: // Tabulated probability tables
          {
            // continuous anisotropic spectra
            TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
            for (int i = 0, eN2 = tab2->GetN2(); i != eN2; ++i) {
              TNudyEndfTab1 *tab = (TNudyEndfTab1 *)recIter.Next();
              GenerateProb(tab);
            }
            if (fNI > 0) sec->RemoveObj(tab2);
          } break;
        } 
      } break;
    } // end switch
  } // end of while loop
} // end loop class

TNudyEndfPhAng::~TNudyEndfPhAng()
{
  fNbt1.shrink_to_fit();
  fInt1.shrink_to_fit();
  fCosFile4.shrink_to_fit();
  fCosPdfFile4.shrink_to_fit();
  fEin.shrink_to_fit();
  fLegendCoef.shrink_to_fit();
}

double TNudyEndfPhAng::RecursionLinearLeg1D(double x1, double x2, double pdf1, double pdf2)
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
double TNudyEndfPhAng::RecursionLinearProb(double x1, double x2, double pdf1, double pdf2)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2) return 0;
  pdf            = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNr, fCosFile4, fCosPdfFile4, fNp, mid);
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
//-----------------------------------------------------------------------------------------
void TNudyEndfPhAng::GenerateProb(TNudyEndfTab1 *tab)
{
  fNr = tab->GetNR();
  fNp = tab->GetNP();
  for (int i = 0; i < fNr; i++) {
    fNbt1.push_back(tab->GetNBT(i));
    fInt1.push_back(tab->GetINT(i));
  }
  for (int j = 0; j < fNp; j++) {
    fCosFile4.push_back(tab->GetX(j));
    fCosPdfFile4.push_back(tab->GetY(j));
  }
  if (fNp > 2) {
    for (int l = 0; l < fNp - 1; l++) {
      RecursionLinearProb(fCosFile4[l], fCosFile4[l + 1], fCosPdfFile4[l], fCosPdfFile4[l + 1]);
    }
  }
  TNudyCore::Instance()->Sort(fCosFile4, fCosPdfFile4);
  double energy = tab->GetC2();
  ModifyTab1(tab, fCosFile4, fCosPdfFile4, energy);
  fCosFile4.clear();
  fCosPdfFile4.clear();
  fNbt1.clear();
  fInt1.clear();  
}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfPhAng::CreateTab2(TNudyEndfTab2 *secTab2, int &NE)
{
  NudyPhysics::TNudyEndfSigma::CreateTab2(secTab2, NE);
}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfPhAng::CreateTab1(TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, int mtf[], int law, double &x)
{
  NudyPhysics::TNudyEndfSigma::CreateTab1(secTab1, x1, x2, mtf, law, x);
}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfPhAng::ModifyTab1(TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, double &x)
{
  NudyPhysics::TNudyEndfSigma::ModifyTab1(secTab1, x1, x2, x);
}  
