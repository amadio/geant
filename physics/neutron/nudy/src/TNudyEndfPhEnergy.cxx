// 	This class is reconstructing probability tables for Energy distribution
// 	of the secondatries
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: March 24, 2016

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
#include "TRandom3.h"
#endif

    TNudyEndfPhEnergy::TNudyEndfPhEnergy()
{
}

//______________________________________________________________________________
TNudyEndfPhEnergy::TNudyEndfPhEnergy(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    TIter recIter(sec->GetRecords());
    for (int k = 0; k < sec->GetN1(); k++) {
      TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
      int MT              = sec->GetMT();
      // int NC              = sec->GetN1();
      fMtNumbers.push_back(MT);
      int LF = tab1->GetL2();
      //std::cout << " LF = " << LF << " MT " << MT << "  NC " << NC << std::endl;
      fNR = tab1->GetN1();
      fNP = tab1->GetN2();
      //****************************************************************************
      // arbitrary tabulated function
      if (LF == 1) {
        for (int cr = 0; cr < fNR; cr++) {
          fNbt1.push_back(tab1->GetNBT(cr));
          fInt1.push_back(tab1->GetINT(cr));
        }
        for (int crs = 0; crs < fNP; crs++) {
          fE1.push_back(tab1->GetX(crs));
          fP1.push_back(tab1->GetY(crs));
        }
        TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
        fNr2                 = tab2->GetN1();
        fNp2                 = tab2->GetN2();

        for (int cr = 0; cr < fNr2; cr++) {
          fNbt2.push_back(tab2->GetNBT(cr));
          fInt2.push_back(tab2->GetINT(cr));
        }
        for (int cr = 0; cr < fNp2; cr++) {
          TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
          fEin.push_back(tab12->GetC2());
          fNr3 = tab12->GetNR();
          fNp3 = tab12->GetNP();
          for (int i = 0; i < fNr3; i++) {
            fNbt3.push_back(tab12->GetNBT(i));
            fInt3.push_back(tab12->GetINT(i));
          }
          for (int crs = 0; crs < fNp3; crs++) {
            fEnergyFile5.push_back(tab12->GetX(crs));
            fEnergyPdfFile5.push_back(tab12->GetY(crs));
// 	    std::cout << " E = " << tab12->GetX(crs) <<"  "<< tab12->GetY(crs) << std::endl ;
          }
          for (int cr = 0; cr < fNp3 - 1; cr++) {
            // std::cout << fEnergyFile5[cr] <<"  "<< fEnergyPdfFile5[cr] << std::endl;
//             RecursionLinearFile5Prob(fEnergyFile5[cr], fEnergyFile5[cr + 1], fEnergyPdfFile5[cr], fEnergyPdfFile5[cr + 1]);
          }
          FillPdf1D();
          fNbt3.clear();
          fInt3.clear();
        }
        fNbt1.clear();
        fInt1.clear();
        fNbt2.clear();
        fInt2.clear();
        fE1.clear();
        //****************************************************************************
        // general evaporation spectrum
      }
    }
    fEin2D.push_back(fEin);
    fFrac2D.push_back(fP1);
    fEne3D.push_back(fEne2D);
    fPdf3D.push_back(fPdf2D);
    fCdf3D.push_back(fCdf2D);
    fEin.clear();
    fEne2D.clear();
    fPdf2D.clear();
    fCdf2D.clear();
    fP1.clear();
  }
  fMt5Values.push_back(fMtNumbers);
  fEnergy5OfMts.push_back(fEin2D);
  fFraction5OfMts.push_back(fFrac2D);
  fEnergyOut5OfMts.push_back(fEne3D);
  fEnergyPdf5OfMts.push_back(fPdf3D);
  fEnergyCdf5OfMts.push_back(fCdf3D);
  fMtNumbers.clear();
  fEin2D.clear();
  fEne3D.clear();
  fPdf3D.clear();
  fCdf3D.clear();
  fFrac2D.clear();
}

TNudyEndfPhEnergy::~TNudyEndfPhEnergy()
{
  fMtNumbers.shrink_to_fit();
  fE1.shrink_to_fit();
  fP1.shrink_to_fit();
  fE2.shrink_to_fit();
  fP2.shrink_to_fit();
  fE3.shrink_to_fit();
  fP3.shrink_to_fit();
  INorm.shrink_to_fit();
  fNbt1.shrink_to_fit();
  fInt1.shrink_to_fit();
  fNbt2.shrink_to_fit();
  fInt2.shrink_to_fit();
  fNbt3.shrink_to_fit();
  fInt3.shrink_to_fit();
  fEnergyFile5.shrink_to_fit();
  fEnergyPdfFile5.shrink_to_fit();
  fEnergyCdfFile5.shrink_to_fit();
  fEin.shrink_to_fit();
  fEneE.shrink_to_fit();
  fCdf.shrink_to_fit();
  fPdf.shrink_to_fit();
  fEne2D.shrink_to_fit();
  fFrac2D.shrink_to_fit();
  fCdf2D.shrink_to_fit();
  fPdf2D.shrink_to_fit();
  fEin2D.shrink_to_fit();
  fEne3D.shrink_to_fit();
  fCdf3D.shrink_to_fit();
  fPdf3D.shrink_to_fit();
}

double TNudyEndfPhEnergy::RecursionLinearFile5Prob(double x1, double x2, double pdf1, double pdf2)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2) return 0;
  pdf            = TNudyCore::Instance()->Interpolate(fNbt3, fInt3, fNr3, fEnergyFile5, fEnergyPdfFile5, fNp3, mid);
  double pdfmid1 = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  if (fabs((pdf - pdfmid1) / pdfmid1) <= 1E-3) {
    return 0;
  }
  fEnergyFile5.push_back(mid);
  fEnergyPdfFile5.push_back(pdf);
  RecursionLinearFile5Prob(x1, mid, pdf1, pdf);
  RecursionLinearFile5Prob(mid, x2, pdf, pdf2);
  return 0;
}

double TNudyEndfPhEnergy::RecursionLinearFile5GenEva(double x1, double x2, double pdf1, double pdf2, double energy)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2) return 0;
  double gx            = 0.0;
  double pe            = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, fE1, fP1, fNP, energy);
  double thetae        = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, fP2, fNp2, energy);
  if (thetae > 0.0) gx = TNudyCore::Instance()->Interpolate(fNbt3, fInt3, fNr3, fE3, fP3, fNp3, mid / thetae);
  pdf                  = 0.0;
  pdf                  = pe * gx;
  double pdfmid1       = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  if (fabs((pdf - pdfmid1) / pdfmid1) <= 1E-3) {
    return 0;
  }
  fEnergyFile5.push_back(mid);
  fEnergyPdfFile5.push_back(pdf);
  RecursionLinearFile5GenEva(x1, mid, pdf1, pdf, energy);
  RecursionLinearFile5GenEva(mid, x2, pdf, pdf2, energy);
  return 0;
}

double TNudyEndfPhEnergy::RecursionLinearFile5Maxwell(double x1, double x2, double pdf1, double pdf2, double energy)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2) return 0;
  double pe        = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, fE1, fP1, fNP, energy);
  double I         = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, INorm, fNp2, energy);
  double thetae    = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, fP2, fNp2, energy);
  pdf              = 0.0;
  if (I > 0.0) pdf = pe * sqrt(mid) * exp(-mid / thetae) / I;
  double pdfmid1   = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  if (fabs((pdf - pdfmid1) / pdfmid1) <= 1E-3) {
    return 0;
  }
  fEnergyFile5.push_back(mid);
  fEnergyPdfFile5.push_back(pdf);
  RecursionLinearFile5Maxwell(x1, mid, pdf1, pdf, energy);
  RecursionLinearFile5Maxwell(mid, x2, pdf, pdf2, energy);
  return 0;
}

double TNudyEndfPhEnergy::RecursionLinearFile5Watt(double x1, double x2, double pdf1, double pdf2, double energy)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2) return 0;
  double pe        = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, fE1, fP1, fNP, energy);
  double I         = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, INorm, fNp2, energy);
  double ae        = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, fP2, fNp2, energy);
  double be        = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE3, fP3, fNp3, energy);
  pdf              = 0.0;
  if (I > 0.0) pdf = pe * sinh(sqrt(be * mid)) * exp(-mid / ae) / I;
  double pdfmid1   = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  if (fabs((pdf - pdfmid1) / pdfmid1) <= 1E-3) {
    return 0;
  }
  fEnergyFile5.push_back(mid);
  fEnergyPdfFile5.push_back(pdf);
  RecursionLinearFile5Watt(x1, mid, pdf1, pdf, energy);
  RecursionLinearFile5Watt(mid, x2, pdf, pdf2, energy);
  return 0;
}
void TNudyEndfPhEnergy::FillPdf1D()
{
  TNudyCore::Instance()->Sort(fEnergyFile5, fEnergyPdfFile5);
  TNudyCore::Instance()->ThinningDuplicate(fEnergyFile5, fEnergyPdfFile5);
  TNudyCore::Instance()->cdfGenerateT(fEnergyFile5, fEnergyPdfFile5, fEnergyCdfFile5);
  for (unsigned long i = 0; i < fEnergyFile5.size(); i++) {
    if (fEnergyPdfFile5[i] > 1E-15) {
      fEneE.push_back(fEnergyFile5[i]);
      fPdf.push_back(fEnergyPdfFile5[i]);
      fCdf.push_back(fEnergyCdfFile5[i]);
    }
  }
  fEne2D.push_back(fEneE);
  fPdf2D.push_back(fPdf);
  fCdf2D.push_back(fCdf);
  fEnergyFile5.clear();
  fEnergyPdfFile5.clear();
  fEnergyCdfFile5.clear();
  fEneE.clear();
  fPdf.clear();
  fCdf.clear();
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfPhEnergy::GetEnergy5(int ielemId, int mt, double energyK)
{
  fRnd  = new TRandom3(0);
  int i = -1;
  // std::cout<<"mt "<< mt <<"  "<< fMt5Values[ielemId].size() <<"  "<<energyK << std::endl;
  for (unsigned int l = 0; l < fMt5Values[ielemId].size(); l++) {
    if (fMt5Values[ielemId][l] == mt) {
      i = l;
      break;
    }
  }
  // std::cout<<"i "<< i <<"  "<< fMt5Values[ielemId][i] <<"  "<<fEnergy5OfMts[ielemId][i].size() << std::endl;
  if (i < 0) return 99;
  int min = 0;
  int max = fEnergy5OfMts[ielemId][i].size() - 1;
  int mid = 0;
  if (energyK <= fEnergy5OfMts[ielemId][i][min])
    min = 0;
  else if (energyK >= fEnergy5OfMts[ielemId][i][max])
    min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < fEnergy5OfMts[ielemId][i][mid])
        max = mid;
      else
        min = mid;
    }
  }
  if (min < 0) min = 0;
  // for(unsigned int kk = 0; kk < fEnergy5OfMts[ielemId][i].size(); kk++)
  // std::cout << min <<"  "<< fEnergy5OfMts[ielemId][i][kk] <<"  "<< fEnergy5OfMts[ielemId][i].size() << std::endl;
  double fraction = (energyK - fEnergy5OfMts[ielemId][i][min]) / 
                    (fEnergy5OfMts[ielemId][i][min + 1] - fEnergy5OfMts[ielemId][i][min]);
  double rnd1              = fRnd->Uniform(1);
  double rnd2              = fRnd->Uniform(1);
  if (rnd2 < fraction) min = min + 1;
  // std::cout<<"min "<< min <<"  "<< fEnergy5OfMts[ielemId][i][min] << std::endl;
  int k    = 0;
  int size = fEnergyCdf5OfMts[ielemId][i][min].size();
  for (unsigned int j = 1; j < fEnergyPdf5OfMts[ielemId][i][min].size(); j++) {
    if (rnd1 <= fEnergyCdf5OfMts[ielemId][i][min][j]) {
      k                    = j - 1;
      if (k >= size - 2) k = size - 2;
      break;
    }
  }
  // for (unsigned int j1 = 0; j1 < fEnergyOut5OfMts[ielemId][i][min].size(); j1++){
  // std::cout<< fEnergyOut5OfMts[ielemId][i][min][j1] <<"  "<< fEnergyPdf5OfMts[ielemId][i][min][j1] <<std::endl;
  //}
  double plk = (fEnergyPdf5OfMts[ielemId][i][min][k + 1] - fEnergyPdf5OfMts[ielemId][i][min][k]) /
               (fEnergyOut5OfMts[ielemId][i][min][k + 1] - fEnergyOut5OfMts[ielemId][i][min][k]);
  double plk2 = fEnergyPdf5OfMts[ielemId][i][min][k] * fEnergyPdf5OfMts[ielemId][i][min][k];

  double edes = 0;
  if (plk != 0)
    edes = 
        fEnergyOut5OfMts[ielemId][i][min][k] +
        (sqrt(plk2 + 2 * plk * (rnd1 - fEnergyCdf5OfMts[ielemId][i][min][k])) - fEnergyPdf5OfMts[ielemId][i][min][k]) /
        plk;
  return edes;
}

//------------------------------------------------------------------------------------------------------
double TNudyEndfPhEnergy::GetDelayedFraction(int ielemId, int mt, double energyK)
{
  fRnd  = new TRandom3(0);
  int i = -1;
  // std::cout<<"mt "<< mt <<"  "<< fMt5Values[ielemId].size() <<"  "<<energyK << std::endl;
  for (unsigned int l = 0; l < fMt5Values[ielemId].size(); l++) {
    if (fMt5Values[ielemId][l] == mt) {
      i = l;
      break;
    }
  }
  // std::cout<<"i "<< i <<"  "<< fMt5Values[ielemId][i] << std::endl;
  if (i < 0) return 99;
  int min = 0;
  int max = fEnergy5OfMts[ielemId][i].size() - 1;
  int mid = 0;
  if (energyK <= fEnergy5OfMts[ielemId][i][min])
    min = 0;
  else if (energyK >= fEnergy5OfMts[ielemId][i][max])
    min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < fEnergy5OfMts[ielemId][i][mid])
        max = mid;
      else
        min = mid;
    }
  }
  return fFraction5OfMts[ielemId][i][min];
}