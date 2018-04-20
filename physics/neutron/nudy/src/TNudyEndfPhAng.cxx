// 	This class is reconstructing probability tables for Angular distribution
// 	of the secondatries
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: March 24, 2016

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
#include "TRandom3.h"
#endif

    TNudyEndfPhAng::TNudyEndfPhAng()
    : TNudyEndfRecoPoint()
{
}
//______________________________________________________________________________
TNudyEndfPhAng::TNudyEndfPhAng(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    TIter recIter(sec->GetRecords());
    // TNudyEndfCont *header = (TNudyEndfCont *)recIter.Next();
    int MT                = sec->GetMT();
    fMtNumbers.push_back(MT);
    int LTT 		 = sec->GetL2();
    int LI  		 = sec->GetL1();
    /*
     * int NK  		 = sec->GetN1();
     * printf("NK = %d LTT = %d LI = %d \n",NK, LTT, LI);
     */
    // Legendre polynomial coefficients
    if (LTT == 1 && LI == 0) {
      TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
      for (int i = 0; i < tab2->GetN2(); i++) {
        TNudyEndfList *tab = (TNudyEndfList *)recIter.Next();
        fEin.push_back(tab->GetC2());
        //std::cout<<"energy "<< tab->GetC2() << std::endl;
        for (int j = 0; j < tab->GetNPL(); j++) {
          fLegendCoef1.push_back(tab->GetLIST(j));
        }
        fLegendCoef.push_back(fLegendCoef1);
        fLegendCoef1.clear();
      }
      for (unsigned long i = 0; i < fEin.size(); i++) {
        // printf("Ein = %e\n", fEin[i]);
        int k1     = 0;
        double fme = 0.0;
        do {
          fme      = 0.5;
          double x = -1. + k1 * 0.02;
          for (unsigned long j = 0; j < fLegendCoef[i].size(); j++) {
            double leg = ROOT::Math::legendre(j + 1, x);
            fme += 0.5 * (2. * (j + 1) + 1.) * fLegendCoef[i][j] * leg;
            // printf("a%d = %e leg= %e\n", j, fLegendCoef[i].At(j),leg);
          }
          if (fme > 0.0) {
            fCosFile4.push_back(x);
            fCosPdfFile4.push_back(fme);
          }
          // printf("%e %e\n", x, fme);
          k1++;
        } while (k1 < 101);
        for (int l = 0; l < 100; l++) {
          RecursionLinearLeg(i, fCosFile4[l], fCosFile4[l + 1], fCosPdfFile4[l], fCosPdfFile4[l + 1]);
        }
        FillPdf1D();
      }
      FillPdf2D();
      fLegendCoef.clear();
      // Tabulated probability tables
    } else if (LTT == 2 && LI == 0) {
      TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
      for (int i = 0; i < tab2->GetN2(); i++) {
        TNudyEndfTab1 *tab = (TNudyEndfTab1 *)recIter.Next();
        fEin.push_back(tab->GetC2());
        fNr = tab->GetNR();
        fNp = tab->GetNP();
        // std::cout<<"energy "<< tab->GetC2() << std::endl;
        for (int il = 0; il < tab->GetNR(); il++) {
          fNbt1.push_back(tab->GetNBT(il));
          fInt1.push_back(tab->GetINT(il));
        }
        for (int j = 0; j < tab->GetNP(); j++) {
          fCosFile4.push_back(tab->GetX(j));
          fCosPdfFile4.push_back(tab->GetY(j));
        }
        for (int l = 0; l < tab->GetNP() - 1; l++) {
          RecursionLinearProb(fCosFile4[l], fCosFile4[l + 1], fCosPdfFile4[l], fCosPdfFile4[l + 1]);
        }
        FillPdf1D();
        fNbt1.clear();
        fInt1.clear();
      }
      FillPdf2D();
      // Low energy given by legendre polynomial and high energy by tabulated probability tables
    } else if (LTT == 3 && LI == 0) {
      TNudyEndfTab2 *lowE = (TNudyEndfTab2 *)recIter.Next();
      for (int i = 0; i < lowE->GetN2(); i++) {
        TNudyEndfList *tab = (TNudyEndfList *)recIter.Next();
        fEin.push_back(tab->GetC2());
        for (int j = 0; j < tab->GetNPL(); j++) {
          fLegendCoef1.push_back(tab->GetLIST(j));
        }
        fLegendCoef.push_back(fLegendCoef1);
        fLegendCoef1.clear();
      }
      for (unsigned long i = 0; i < fEin.size(); i++) {
        // printf("Ein = %e\n", fEin[i]);
        int k1     = 0;
        double fme = 0.0;
        do {
          fme      = 0.5;
          double x = -1. + k1 * 0.02;
          for (unsigned long j = 0; j < fLegendCoef[i].size(); j++) {
            double leg = ROOT::Math::legendre(j + 1, x);
            fme += 0.5 * (2. * (j + 1) + 1.) * fLegendCoef[i][j] * leg;
            //printf("a%e = %e leg= %e\n", x, fLegendCoef[i][j],leg);
          }
          if (fme > 0.0) {
            fCosFile4.push_back(x);
            fCosPdfFile4.push_back(fme);
          }
          // printf("%e %e\n", x, fme);
          k1++;
        } while (k1 < 101);

        for (int l = 0; l < 100; l++) {
          RecursionLinearLeg(i, fCosFile4[l], fCosFile4[l + 1], fCosPdfFile4[l], fCosPdfFile4[l + 1]);
        }
        FillPdf1D();
      }
      fLegendCoef.clear();
      TNudyEndfTab2 *highE = (TNudyEndfTab2 *)recIter.Next();
      for (int i = 0; i < highE->GetN2(); i++) {
        TNudyEndfTab1 *tab = (TNudyEndfTab1 *)recIter.Next();
        fEin.push_back(tab->GetC2());
        // std::cout <<"energy "<< fEin[fEin.size() - 1] << std::endl;
        fNr = tab->GetNR();
        fNp = tab->GetNP();
        for (int i = 0; i < tab->GetNR(); i++) {
          fNbt1.push_back(tab->GetNBT(i));
          fInt1.push_back(tab->GetINT(i));
        }
        for (int j = 0; j < tab->GetNP(); j++) {
          fCosFile4.push_back(tab->GetX(j));
          fCosPdfFile4.push_back(tab->GetY(j));
        }
        if (fNp > 2) {
          for (int l = 0; l < tab->GetNP() - 1; l++) {
            RecursionLinearProb(fCosFile4[l], fCosFile4[l + 1], fCosPdfFile4[l], fCosPdfFile4[l + 1]);
          }
        }
        FillPdf1D();
        fNbt1.clear();
        fInt1.clear();
      }
      FillPdf2D();
    } else if (LTT == 0 && LI == 1) {
      fEin.push_back(1E-14);
      fEin.push_back(1.5E8);
      for (int j = 0; j < 2; j++) {
        fCosFile4.push_back(1);
        fCosPdfFile4.push_back(0.5);
        fCosFile4.push_back(-1);
        fCosPdfFile4.push_back(0.5);
        FillPdf1D();
      }
      FillPdf2D();
    }
    // Low energy given by legendre polynomial and high energy by tabulated probability tables
  } // end while loop
  fMt4Values.push_back(fMtNumbers);
  fMt4Lct.push_back(fMtLct);
  fEnergy4OfMts.push_back(fEin2D);
  fCos4OfMts.push_back(fCos3D);
  fCosPdf4OfMts.push_back(fPdf3D);
  fCosCdf4OfMts.push_back(fCdf3D);
  std::vector<int>().swap(fMtNumbers);
  std::vector<int>().swap(fMtLct);
  fEin2D.clear();
  fCos3D.clear();
  fPdf3D.clear();
  fCdf3D.clear();
  /*
  std::cout<<"energies "<< fEnergy4OfMts.size() << std::endl;
  for(unsigned long i = 0; i < fEnergy4OfMts.size(); i++){
    std::cout <<" mt "<<fMt4Values[0][i]<<" size "<< fEnergy4OfMts[i].size() << std::endl;
    for(unsigned long j = 0; j < fEnergy4OfMts[i].size(); j++){
      std::cout << fEnergy4OfMts[i][j] << std::endl;
      for(unsigned long k =0; k < fCosPdf4OfMts[i][j].size()/2; k++){
  std::cout << fCosPdf4OfMts[i][j][2*k] <<"  "<< fCosPdf4OfMts[i][j][2*k + 1] <<"  "<< fCosCdf4OfMts[i][j][2*k + 1]<<
  std::endl;
      }
    }
  }
  */
} // end class

TNudyEndfPhAng::~TNudyEndfPhAng()
{
  fMtLct.shrink_to_fit();
  fMtNumbers.shrink_to_fit();
  fNbt1.shrink_to_fit();
  fInt1.shrink_to_fit();
  fCosFile4.shrink_to_fit();
  fCosPdfFile4.shrink_to_fit();
  fCosCdfFile4.shrink_to_fit();
  fEin.shrink_to_fit();
  fCos4.shrink_to_fit();
  fCdf.shrink_to_fit();
  fPdf.shrink_to_fit();
  fLegendCoef1.shrink_to_fit();
  fCos2D.shrink_to_fit();
  fLegendCoef.shrink_to_fit();
  fCdf2D.shrink_to_fit();
  fPdf2D.shrink_to_fit();
  fEin2D.shrink_to_fit();
  fCos3D.shrink_to_fit();
  fCdf3D.shrink_to_fit();
  fPdf3D.shrink_to_fit();
}

double TNudyEndfPhAng::RecursionLinearLeg(int i, double x1, double x2, double pdf1, double pdf2)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2) return 0;
  //  std::cout <<" beg   "<< x1 <<"  "<< x2 <<"  "<< pdf1 <<"  "<< pdf2 << std::endl;
  for (unsigned long j = 0; j < fLegendCoef[i].size(); j++) {
    double leg = ROOT::Math::legendre(j + 1, mid);
    pdf += 0.5 * (2. * (j + 1) + 1.) * fLegendCoef[i][j] * leg;
  }
  double pdfmid1 = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  // std::cout <<x1<<"  " <<x2<<"  "<< pdf1 <<"  "<< pdf2 << std::endl;
  // std::cout <<mid<<"  "<< fabs((pdf - pdfmid1)/pdf) <<"  "<< pdf <<"  "<< pdfmid1 << std::endl;
  if (pdf > 0 && fabs((pdf - pdfmid1) / pdf) <= 2E-4) {
    return 0;
  }
  fCosFile4.push_back(mid);
  fCosPdfFile4.push_back(pdf);
  RecursionLinearLeg(i, x1, mid, pdf1, pdf);
  RecursionLinearLeg(i, mid, x2, pdf, pdf2);
  return 0;
}
//--------------------------------------------------------------------------------------
double TNudyEndfPhAng::RecursionLinearProb(double x1, double x2, double pdf1, double pdf2)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2) return 0;
  // std::cout <<" beg   "<< x1 <<"  "<< x2 <<"  "<< pdf1 <<"  "<< pdf2 << std::endl;
  pdf            = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNr, fCosFile4, fCosPdfFile4, fNp, mid);
  double pdfmid1 = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  //  std::cout << mid <<" linear "<< pdf <<"  "<< pdfmid1 << std::endl;

  if (pdf > 0 && fabs((pdf - pdfmid1) / pdf) <= 2E-4) {
    return 0;
  }
  fCosFile4.push_back(mid);
  fCosPdfFile4.push_back(pdf);
  RecursionLinearProb(x1, mid, pdf1, pdf);
  RecursionLinearProb(mid, x2, pdf, pdf2);
  return 0;
}
//-----------------------------------------------------------------------------------------
void TNudyEndfPhAng::FillPdf1D()
{
  TNudyCore::Instance()->Sort(fCosFile4, fCosPdfFile4);
  TNudyCore::Instance()->cdfGenerateT(fCosFile4, fCosPdfFile4, fCosCdfFile4);
  for (unsigned long i = 0; i < fCosFile4.size(); i++) {
    fCos4.push_back(fCosFile4[i]);
    fPdf.push_back(fCosPdfFile4[i]);
    fCdf.push_back(fCosCdfFile4[i]);
    // std::cout<<fCosFile4[i] <<"  "<<fCosPdfFile4[i] <<"  "<< fCosCdfFile4[i] << std::endl;
  }
  fCos2D.push_back(fCos4);
  fPdf2D.push_back(fPdf);
  fCdf2D.push_back(fCdf);
  fCosFile4.clear();
  fCosPdfFile4.clear();
  fCosCdfFile4.clear();
  fCos4.clear();
  fPdf.clear();
  fCdf.clear();
}
//--------------------------------------------------------------------------------------
void TNudyEndfPhAng::FillPdf2D()
{
  fEin2D.push_back(fEin);
  fCos3D.push_back(fCos2D);
  fPdf3D.push_back(fPdf2D);
  fCdf3D.push_back(fCdf2D);
  fEin.clear();
  fCos2D.clear();
  fPdf2D.clear();
  fCdf2D.clear();
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfPhAng::GetCos4(int ielemId, int mt, double energyK)
{
  fRnd  = new TRandom3(0);
  int i = -1;
  for (unsigned int l = 0; l < fMt4Values[ielemId].size(); l++) {
    if (fMt4Values[ielemId][l] == mt) {
      i = l;
      break;
    }
  }
  // std::cout<<"i "<< i<< "  " << fEnergy4OfMts[ielemId][i].size() << std::endl;
  int min = 0;
  int max = fEnergy4OfMts[ielemId][i].size() - 1;
  int mid = 0;
  if (energyK <= fEnergy4OfMts[ielemId][i][min])
    min = 0;
  else if (energyK >= fEnergy4OfMts[ielemId][i][max])
    min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < fEnergy4OfMts[ielemId][i][mid])
        max = mid;
      else
        min = mid;
    }
  }
  // std::cout<<" min "<< min << std::endl;
  double fraction = (energyK - fEnergy4OfMts[ielemId][i][min]) / 
                    (fEnergy4OfMts[ielemId][i][min + 1] - fEnergy4OfMts[ielemId][i][min]);
  // std::cout <<" fraction "<< fraction <<"  "<< energyK <<"  "<< fEnergy4OfMts[ielemId][i][min] << std::endl;
  double rnd1              = fRnd->Uniform(1);
  double rnd2              = fRnd->Uniform(1);
  if (rnd2 < fraction) min = min + 1;
  int k                    = 0;
  int size = fCosCdf4OfMts[ielemId][i][min].size();
  for (int j = 1; j < size; j++) {
    if (rnd1 <= fCosCdf4OfMts[ielemId][i][min][j]) {
      k                    = j - 1;
      if (k >= size - 2) k = size - 2;
      break;
    }
  }
  double plk = (fCosPdf4OfMts[ielemId][i][min][k + 1] - fCosPdf4OfMts[ielemId][i][min][k]) /
               (fCos4OfMts[ielemId][i][min][k + 1] - fCos4OfMts[ielemId][i][min][k]);
  double plk2       = fCosPdf4OfMts[ielemId][i][min][k] * fCosPdf4OfMts[ielemId][i][min][k];
  double plsq       = plk2 + 2 * plk * (rnd1 - fCosCdf4OfMts[ielemId][i][min][k]);
  double Ang        = 0;
  if (plk == 0) Ang = fCos4OfMts[ielemId][i][min][k];
  if (plk != 0 && plsq > 0) {
    Ang = fCos4OfMts[ielemId][i][min][k] + (sqrt(std::fabs(plsq)) - fCosPdf4OfMts[ielemId][i][min][k]) / plk;
  } else {
    // std::cout <<"uniform angle "<< mt << std::endl;
    Ang = 2 * rnd1 - 1;
  }
  return Ang;
}
//________________________________________________________________________________________________________________
int TNudyEndfPhAng::GetCos4Lct(int ielemId, int mt)
{
  int i = 0;
  for (unsigned int l = 0; l < fMt4Values[ielemId].size(); l++) {
    if (fMt4Values[ielemId][l] == mt) {
      i = l;
      break;
    }
  }
  return fMt4Lct[ielemId][i];
}