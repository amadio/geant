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
#include "Geant/TNudyEndfEnergy.h"
#include "Math/SpecFuncMathMore.h"
#include "TMath.h"

using namespace Nudy;
using namespace NudyPhysics;

#ifdef USE_ROOT
ClassImp(TNudyEndfEnergy)
#endif

TNudyEndfEnergy::TNudyEndfEnergy() : TNudyEndfSigma()
{
}
//______________________________________________________________________________
TNudyEndfEnergy::TNudyEndfEnergy(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    TIter recIter(sec->GetRecords());
    int NK     = sec->GetN1();
    for (int k = 0; k < NK; k++) {
      TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
      int MT              = sec->GetMT();
      fMtNumbers.push_back(MT);
      int LF = tab1->GetL2();
      mtf[0] = sec->GetMAT();
      mtf[1] = sec->GetMT();
      mtf[2] = sec->GetMF();
      
      switch (LF){ 
        case 1: // arbitrary tabulated functions
        {
          ProcessTab1(tab1, fNR, fNP, fNbt1, fInt1, fE1, fP1);
          TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
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
          }break;
        case 5: // general evaporation spectrum
        {  
          double u = tab1->GetC1();
          ProcessTab1(tab1, fNR, fNP, fNbt1, fInt1, fE1, fP1);
          tab1->SetCont(tab1->GetC1(),tab1->GetC2(),tab1->GetL1(),1,fNR,fNP);
          for (int i = 0; i < fNR; i++) {
            tab1->SetNBT(fNbt1[i],i);
            tab1->SetINT(fInt1[i],i);
          }
          for (int crs = 0; crs < fNP; crs++) {
            tab1->SetX(fE1[crs],crs);
            tab1->SetY(fP1[crs],crs);
          }
          TNudyEndfTab1 *tab11 = (TNudyEndfTab1 *)recIter.Next();
          ProcessTab1(tab11, fNr2, fNp2, fNbt2, fInt2, fE2, fP2);
          // creating Tab2 parameters and inserting after tab1 as required by LF = 1
          TNudyEndfTab2 *secTab2 = new TNudyEndfTab2();
          CreateTab2(secTab2,fNp2);
          sec->AddAfter(tab1,secTab2);
          sec->RemoveObj(tab11);
          
          TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
          ProcessTab1(tab12, fNr3, fNp3, fNbt3, fInt3, fE3, fP3);
          sec->RemoveObj(tab12);
          TNudyEndfTab1 *secTab1(new TNudyEndfTab1[fNp2]);
          double energy, eout;
          for (int i = 0; i < fNp2; i++) {
            energy             = fE2[i];
            if (i == 0)  energy = energy * 1.001;
            eout               = 1E-05;
            double EHLimit     = energy - u;
            if (mtf[1] == 18 || mtf[1] == 19 || mtf[1] == 20 || mtf[1] == 21 || mtf[1] == 38 || mtf[1] == 455) {
              EHLimit = fabs(u); }
            int    FlagElimit  = -1 ;
            double sumprob     = 0.0;
            do {
              if (eout >= EHLimit - 1E-5) FlagElimit = 1;
              double gx     = 0.0;
              double pe     = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, fE1, fP1, fNP, energy);
              double thetae = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, fP2, fNp2, energy);
              if (thetae > 0.0)
                gx       = TNudyCore::Instance()->Interpolate(fNbt3, fInt3, fNr3, fE3, fP3, fNp3, eout / thetae);
              double ppe = 0.0;
              ppe        = pe * gx;
              if (ppe >= 0.0) {
                fEnergyFile5.push_back(eout);
                fEnergyPdfFile5.push_back(ppe);
                sumprob += ppe;
              }
              eout *= 2;
              if (eout >= EHLimit && FlagElimit == -1) {
                eout = EHLimit - 1E-5 ;
              }
            } while (eout < EHLimit);
            if (sumprob >= 0.0) {
              fEin.push_back(energy);
              int size = fEnergyFile5.size();
              for (int cr = 0; cr < size - 1; cr++) {
                RecursionLinearFile5GenEva(fEnergyFile5[cr], fEnergyFile5[cr + 1], fEnergyPdfFile5[cr],
                                          fEnergyPdfFile5[cr + 1], energy);
              }
              TNudyCore::Instance()->Sort(fEnergyFile5, fEnergyPdfFile5);
              // creating Tab1 parameters and inserting after tab1 as required by LF = 1
              CreateTab1(&secTab1[i], fEnergyFile5, fEnergyPdfFile5, mtf, 2, energy);
              if (i == 0) {
                sec->AddAfter(secTab2, &secTab1[i]);
              } else {
                sec->AddAfter(&secTab1[i-1], &secTab1[i]);
              }
              fEnergyFile5.clear();
              fEnergyPdfFile5.clear();
            }
          }
          fNbt1.clear();
          fInt1.clear();
          fE1.clear();
          fP1.clear();
          fNbt2.clear();
          fInt2.clear();
          fE2.clear();
          fP2.clear();
          fNbt3.clear();
          fInt3.clear();
          fE3.clear();
          fP3.clear();
        }break;
        case 7: // simple maxwellian fission spectrum
        { 
          double u = tab1->GetC1();
          ProcessTab1(tab1, fNR, fNP, fNbt1, fInt1, fE1, fP1);
          tab1->SetCont(tab1->GetC1(),tab1->GetC2(),tab1->GetL1(),1,fNR,fNP);
          for (int i = 0; i < fNR; i++) {
            tab1->SetNBT(fNbt1[i],i);
            tab1->SetINT(fInt1[i],i);
          }
          for (int crs = 0; crs < fNP; crs++) {
            tab1->SetX(fE1[crs],crs);
            tab1->SetY(fP1[crs],crs);
          }
          TNudyEndfTab1 *tab11 = (TNudyEndfTab1 *)recIter.Next();
          ProcessTab1(tab11, fNr2, fNp2, fNbt2, fInt2, fE2, fP2);
          // creating Tab2 parameters and inserting after tab1 as required by LF = 1
          TNudyEndfTab2 *secTab2 = new TNudyEndfTab2();
          CreateTab2(secTab2,fNp2);
          sec->AddAfter(tab1,secTab2);
          sec->RemoveObj(tab11);
          for (int crs = 0; crs < fNp2; crs++) {
            double sq  = (fE2[crs] - u) / fP2[crs];
            double sqt = sqrt(sq);
            INorm.push_back(pow(fP2[crs], 1.5) * (0.5 * 1.7724538529055 * erf(sqt) - sqt * exp(-sq)));
          }
          TNudyEndfTab1 *secTab1(new TNudyEndfTab1[fNp2]);
          double energy, eout;
          for (int i = 0; i < fNp2; i++) {
            energy             = fE2[i];
            if (i == 0)  energy = energy * 1.001;
            eout               = 1E-05;
            double EHLimit     = energy - u;
            if (mtf[1] == 18 || mtf[1] == 19 || mtf[1] == 20 || mtf[1] == 21 || mtf[1] == 38 || mtf[1] == 455) {
              EHLimit = fabs(u); }
            double sumprob     = 0.0;
            int    FlagElimit  = -1 ;
            do {
              if (eout >= EHLimit - 1E-5) FlagElimit = 1;
              double pe        = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, fE1, fP1, fNP, energy);
              double I         = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, INorm, fNp2, energy);
              double thetae    = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, fP2, fNp2, energy);
              double ppe       = 0.0;
              if (I > 0.0) ppe = pe * sqrt(eout) * exp(-eout / thetae) / I;
              if (ppe >= 0.0) {
                fEnergyFile5.push_back(eout);
                fEnergyPdfFile5.push_back(ppe);
                sumprob += ppe;
              }
              eout *= 2;
              if (eout >= EHLimit && FlagElimit == -1) {
                eout = EHLimit - 1E-5;
              }
            } while (eout <= EHLimit);
            if (sumprob >= 0.0) {
              fEin.push_back(energy);
              int size = fEnergyFile5.size();
              for (int cr = 0; cr < size - 1; cr++) {
                RecursionLinearFile5Maxwell(fEnergyFile5[cr], fEnergyFile5[cr + 1], fEnergyPdfFile5[cr],
                                            fEnergyPdfFile5[cr + 1], energy);
              }
              TNudyCore::Instance()->Sort(fEnergyFile5, fEnergyPdfFile5);
              // creating Tab1 parameters and inserting after tab1 as required by LF = 1
              CreateTab1(&secTab1[i], fEnergyFile5, fEnergyPdfFile5, mtf, 2, energy);
              if (i == 0) {
                sec->AddAfter(secTab2, &secTab1[i]);
              } else {
                sec->AddAfter(&secTab1[i-1], &secTab1[i]);
              }
              fEnergyFile5.clear();
              fEnergyPdfFile5.clear();
            }
          }
          fNbt1.clear();
          fInt1.clear();
          fE1.clear();
          fP1.clear();
          fNbt2.clear();
          fInt2.clear();
          fE2.clear();
          fP2.clear();
          INorm.clear();
        }break;
        case 9: // evaporation spectrum
        {
          double u = tab1->GetC1();
          ProcessTab1(tab1, fNR, fNP, fNbt1, fInt1, fE1, fP1);
          tab1->SetCont(tab1->GetC1(),tab1->GetC2(),tab1->GetL1(),1,fNR,fNP);
          for (int i = 0; i < fNR; i++) {
            tab1->SetNBT(fNbt1[i],i);
            tab1->SetINT(fInt1[i],i);
          }
          for (int crs = 0; crs < fNP; crs++) {
            tab1->SetX(fE1[crs],crs);
            tab1->SetY(fP1[crs],crs);
          }
          TNudyEndfTab1 *tab11 = (TNudyEndfTab1 *)recIter.Next();
          ProcessTab1(tab11, fNr2, fNp2, fNbt2, fInt2, fE2, fP2);
          // creating Tab2 parameters and inserting after tab1 as required by LF = 1
          TNudyEndfTab2 *secTab2 = new TNudyEndfTab2();
          CreateTab2(secTab2,fNp2);
          sec->AddAfter(tab1,secTab2);
          sec->RemoveObj(tab11);
          for (int crs = 0; crs < fNp2; crs++) {
            INorm.push_back(fP2[crs] * fP2[crs] *
                            (1. - exp(-(fE2[crs] - u) / fP2[crs]) * (1.0 + (fE2[crs] - u) / fP2[crs])));
          }
          TNudyEndfTab1 *secTab1(new TNudyEndfTab1[fNp2]);
          double energy, eout;
          for (int i = 0; i < fNp2; i++) {
            energy             = fE2[i];
            if (i == 0)  energy = energy * 1.001;
            eout               = 1E-05;
            double EHLimit     = energy - u;
            if (mtf[1] == 18 || mtf[1] == 19 || mtf[1] == 20 || mtf[1] == 21 || mtf[1] == 38 || mtf[1] == 455) {
              EHLimit = fabs(u); 
            }
            double sumprob     = 0.0;
            int    FlagElimit  = -1 ;
            do {
              if (eout >= EHLimit - 1E-5) FlagElimit = 1;
              double pe        = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, fE1, fP1, fNP, energy);
              double I         = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, INorm, fNp2, energy);
              double thetae    = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, fP2, fNp2, energy);
              double ppe       = 0.0;
              if (I > 0.0) ppe = pe * eout * exp(-eout / thetae) / I;
              if (ppe >= 0.0) {
                fEnergyFile5.push_back(eout);
                fEnergyPdfFile5.push_back(ppe);
                sumprob += ppe;
              }
              eout *= 2;
              if (eout >= EHLimit && FlagElimit == -1) {
                eout = EHLimit - 1E-5 ;
              }
            } while (eout <= EHLimit);
            if (sumprob >= 0.0) {
              fEin.push_back(energy);
              int size = fEnergyFile5.size();
              for (int cr = 0; cr < size - 1; cr++) {
                RecursionLinearFile5Maxwell(fEnergyFile5[cr], fEnergyFile5[cr + 1], fEnergyPdfFile5[cr],
                                            fEnergyPdfFile5[cr + 1], energy);
              }
              TNudyCore::Instance()->Sort(fEnergyFile5, fEnergyPdfFile5);
              // creating Tab1 parameters and inserting after tab1 as required by LF = 1
              CreateTab1(&secTab1[i], fEnergyFile5, fEnergyPdfFile5, mtf, 2, energy);
              if (i == 0) {
                sec->AddAfter(secTab2, &secTab1[i]);
              } else {
                sec->AddAfter(&secTab1[i-1], &secTab1[i]);
              }
              fEnergyFile5.clear();
              fEnergyPdfFile5.clear();
            }
          }
          fNbt1.clear();
          fInt1.clear();
          fE1.clear();
          fP1.clear();
          fNbt2.clear();
          fInt2.clear();
          fE2.clear();
          fP2.clear();
          INorm.clear();
        }break;
        case 11: // energy dependent watt spectrum
        {
          double u = tab1->GetC1();
          ProcessTab1(tab1, fNR, fNP, fNbt1, fInt1, fE1, fP1);
          tab1->SetCont(tab1->GetC1(),tab1->GetC2(),tab1->GetL1(),1,fNR,fNP);
          for (int i = 0; i < fNR; i++) {
            tab1->SetNBT(fNbt1[i],i);
            tab1->SetINT(fInt1[i],i);
          }
          for (int crs = 0; crs < fNP; crs++) {
            tab1->SetX(fE1[crs],crs);
            tab1->SetY(fP1[crs],crs);
          }
          TNudyEndfTab1 *tab11 = (TNudyEndfTab1 *)recIter.Next();
          ProcessTab1(tab11, fNr2, fNp2, fNbt2, fInt2, fE2, fP2);
          // creating Tab2 parameters and inserting after tab1 as required by LF = 1
          TNudyEndfTab2 *secTab2 = new TNudyEndfTab2();
          CreateTab2(secTab2,fNp2);
          sec->AddAfter(tab1,secTab2);
          sec->RemoveObj(tab11);
          
          TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
          ProcessTab1(tab12, fNr3, fNp3, fNbt3, fInt3, fE3, fP3);
          sec->RemoveObj(tab12);
          for (int crs = 0; crs < fNp3; crs++) {
            double a   = fP2[crs];
            double b   = fP3[crs];
            double eua = (fE3[crs] - u) / a;
            INorm.push_back(0.5 * sqrt(0.25 * PI * a * a * a * b) * exp(0.25 * a * b) *
                            (erf(sqrt(eua) - sqrt(0.25 * a * b)) + erf(sqrt(eua) + sqrt(0.25 * a * b))) -
                            a * exp(-eua) * sinh(sqrt(b * (fE3[crs] - u))));
          }
          TNudyEndfTab1 *secTab1(new TNudyEndfTab1[fNp2]);
          double energy, eout;
          for (int i = 0; i < fNp2; i++) {
            energy             = fE2[i];
            if (i == 0)  energy = energy * 1.001;
            eout               = 1E-05;
            double EHLimit     = energy - u;
            if (mtf[1] == 18 || mtf[1] == 19 || mtf[1] == 20 || mtf[1] == 21 || mtf[1] == 38 || mtf[1] == 455) {
              EHLimit = fabs(u); }
            double sumprob     = 0.0;
            int    FlagElimit  = -1 ;
            do {
              if (eout >= EHLimit - 1E-5) FlagElimit = 1;
              double pe        = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, fE1, fP1, fNP, energy);
              double I         = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, INorm, fNp2, energy);
              double ae        = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, fP2, fNp2, energy);
              double be        = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE3, fP3, fNp3, energy);
              double ppe       = 0.0;
              if (I > 0.0) ppe = pe * sinh(sqrt(be * eout)) * exp(-eout / ae) / I;
              if (ppe >= 0.0) {
                fEnergyFile5.push_back(eout);
                fEnergyPdfFile5.push_back(ppe);
                sumprob += ppe;
              }
              eout *= 2;
              if (eout >= EHLimit && FlagElimit == -1) {
                eout = EHLimit - 1E-5 ;
              }
            } while (eout <= EHLimit);
            if (sumprob >= 0.0) {
              fEin.push_back(energy);
              int size = fEnergyFile5.size();
              for (int cr = 0; cr < size - 1; cr++) {
                RecursionLinearFile5Watt(fEnergyFile5[cr], fEnergyFile5[cr + 1], fEnergyPdfFile5[cr],
                                        fEnergyPdfFile5[cr + 1], energy);
              }
              TNudyCore::Instance()->Sort(fEnergyFile5, fEnergyPdfFile5);
              // creating Tab1 parameters and inserting after tab1 as required by LF = 1
              CreateTab1(&secTab1[i], fEnergyFile5, fEnergyPdfFile5, mtf, 2, energy);
              if (i == 0) {
                sec->AddAfter(secTab2, &secTab1[i]);
              } else {
                sec->AddAfter(&secTab1[i-1], &secTab1[i]);
              }
              fEnergyFile5.clear();
              fEnergyPdfFile5.clear();
            }
          }
          fNbt1.clear();
          fInt1.clear();
          fE1.clear();
          fP1.clear();
          fNbt2.clear();
          fInt2.clear();
          fE2.clear();
          fP2.clear();
          fNbt3.clear();
          fInt3.clear();
          fE3.clear();
          fP3.clear();
          INorm.clear();
        }break;
        case 12: // this law is not tested
        {
          ProcessTab1(tab1, fNR, fNP, fNbt1, fInt1, fE1, fP1);
          tab1->SetCont(tab1->GetC1(),tab1->GetC2(),tab1->GetL1(),1,fNR,fNP);
          for (int i = 0; i < fNR; i++) {
            tab1->SetNBT(fNbt1[i],i);
            tab1->SetINT(fInt1[i],i);
          }
          for (int crs = 0; crs < fNP; crs++) {
            tab1->SetX(fE1[crs],crs);
            tab1->SetY(fP1[crs],crs);
          }
          TNudyEndfTab1 *tab11 = (TNudyEndfTab1 *)recIter.Next();
          ProcessTab1(tab11, fNr2, fNp2, fNbt2, fInt2, fE2, fP2);
          // creating Tab2 parameters and inserting after tab1 as required by LF = 1
          TNudyEndfTab2 *secTab2 = new TNudyEndfTab2();
          CreateTab2(secTab2,fNp2);
          sec->AddAfter(tab1,secTab2);
          sec->RemoveObj(tab11);
          double efl           = tab11->GetC1();
          double efh           = tab11->GetC2();
          for (int crs = 0; crs < fNp2; crs++) {
            fEin.push_back(tab11->GetX(crs));
          }
          TNudyEndfTab1 *secTab1(new TNudyEndfTab1[fNp2]);
          double energy, eout;
          for (int i = 0; i < fNp2; i++) {
            energy             = fE2[i];
            if (i == 0)  energy = energy * 1.001;
            eout                   = 1E-05;
            do {
              double pe  = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, fE1, fP1, fNP, energy);
              double tm  = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, fP2, fNp2, energy);
              double u1l = (sqrt(eout) - sqrt(efl)) * (sqrt(eout) - sqrt(efl)) / tm;
              double u2l = (sqrt(eout) + sqrt(efl)) * (sqrt(eout) + sqrt(efl)) / tm;
              double u1h = (sqrt(eout) - sqrt(efh)) * (sqrt(eout) - sqrt(efh)) / tm;
              double u2h = (sqrt(eout) + sqrt(efh)) * (sqrt(eout) + sqrt(efh)) / tm;
              double e1ul = ROOT::Math::expint(u1l);
              double e2ul = ROOT::Math::expint(u2l);
              double e1uh = ROOT::Math::expint(u1h);
              double e2uh = ROOT::Math::expint(u2h);
              double a1      = 1.5;
              double gamau1l = TMath::Gamma(a1, u1l);
              double gamau2l = TMath::Gamma(a1, u2l);
              double gamau1h = TMath::Gamma(a1, u1h);
              double gamau2h = TMath::Gamma(a1, u2h);
              double gl = (1. / (3 * sqrt(efl * tm))) * (pow(u2l, 1.5) * e2ul - pow(u1l, 1.5) * e1ul + gamau2l - gamau1l);
              double gh = (1. / (3 * sqrt(efh * tm))) * (pow(u2h, 1.5) * e2uh - pow(u1h, 1.5) * e1uh + gamau2h - gamau1h);
              double ppe = 0.5 * pe * (gl + gh);
              fEnergyFile5.push_back(eout);
              fEnergyPdfFile5.push_back(ppe);
              eout *= 2;
            } while (eout < fE1[fNP - 1]);
            // creating Tab1 parameters and inserting after tab1 as required by LF = 1
            CreateTab1(&secTab1[i], fEnergyFile5, fEnergyPdfFile5, mtf, 2, energy);
            if (i == 0) {
              sec->AddAfter(secTab2, &secTab1[i]);
            } else {
              sec->AddAfter(&secTab1[i-1], &secTab1[i]);
            }
            fEnergyFile5.clear();
            fEnergyPdfFile5.clear();
          }
          fNbt1.clear();
          fInt1.clear();
          fE1.clear();
          fP1.clear();
          fNbt2.clear();
          fInt2.clear();
          fE2.clear();
          fP2.clear();
        }
      }
    }
    fEin.clear();
  }
}

TNudyEndfEnergy::~TNudyEndfEnergy()
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
  fEin.shrink_to_fit();
}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfEnergy::ProcessTab2(Nudy::TNudyEndfTab2 *tab2, int &NR,int &NP,rowint &fnbt, rowint &fint)
{
  NudyPhysics::TNudyEndfSigma::ProcessTab2(tab2, NR, NP, fnbt, fint);
}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfEnergy::ProcessTab1(Nudy::TNudyEndfTab1 *tab1,int &NR,int &NP,rowint &fnbt, 
                                  rowint &fint,rowd &x1, rowd &x2)
{  
  NudyPhysics::TNudyEndfSigma::ProcessTab1(tab1, NR, NP, fnbt, fint, x1, x2);
}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfEnergy::ModifyTab1(TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, double &x)
{
  NudyPhysics::TNudyEndfSigma::ModifyTab1(secTab1, x1, x2, x);
}  
// -------------------------------------------------------------------------------------------------------
void TNudyEndfEnergy::CreateTab2(TNudyEndfTab2 *secTab2, int &NE)
{
  NudyPhysics::TNudyEndfSigma::CreateTab2(secTab2, NE);
}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfEnergy::CreateTab1(TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, int mtf[], int law, double &x)
{
  NudyPhysics::TNudyEndfSigma::CreateTab1(secTab1, x1, x2, mtf, law, x);
}
// -------------------------------------------------------------------------------------------------------
double TNudyEndfEnergy::RecursionLinearFile5Prob(double x1, double x2, double pdf1, double pdf2)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2 || pdf2 == 0) return 0;
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
double TNudyEndfEnergy::RecursionLinearFile5GenEva(double x1, double x2, double pdf1, double pdf2, double energy)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2 || pdf2 == 0) return 0;
  double gx            = 0.0;
  double pe            = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, fE1, fP1, fNP, energy);
  double thetae        = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, fP2, fNp2, energy);
  if (thetae > 0.0) gx = TNudyCore::Instance()->Interpolate(fNbt3, fInt3, fNr3, fE3, fP3, fNp3, mid / thetae);
  pdf                  = 0.0;
  pdf                  = pe * gx;
  double pdfmid1       = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  if (fabs((pdf - pdfmid1) / pdfmid1) <= 5E-3) {
    return 0;
  }
  fEnergyFile5.push_back(mid);
  fEnergyPdfFile5.push_back(pdf);
  RecursionLinearFile5GenEva(x1, mid, pdf1, pdf, energy);
  RecursionLinearFile5GenEva(mid, x2, pdf, pdf2, energy);
  return 0;
}
// -------------------------------------------------------------------------------------------------------
double TNudyEndfEnergy::RecursionLinearFile5Maxwell(double x1, double x2, double pdf1, double pdf2, double energy)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2 || pdf2 == 0) return 0;
  double pe        = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, fE1, fP1, fNP, energy);
  double I         = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, INorm, fNp2, energy);
  double thetae    = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, fP2, fNp2, energy);
  pdf              = 0.0;
  if (I > 0.0) pdf = pe * sqrt(mid) * exp(-mid / thetae) / I;
  double pdfmid1   = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  if (fabs((pdf - pdfmid1) / pdfmid1) <= 5E-3) {
    return 0;
  }
  fEnergyFile5.push_back(mid);
  fEnergyPdfFile5.push_back(pdf);
  RecursionLinearFile5Maxwell(x1, mid, pdf1, pdf, energy);
  RecursionLinearFile5Maxwell(mid, x2, pdf, pdf2, energy);
  return 0;
}
// -------------------------------------------------------------------------------------------------------
double TNudyEndfEnergy::RecursionLinearFile5Watt(double x1, double x2, double pdf1, double pdf2, double energy)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2 || pdf2 == 0) return 0;
  double pe        = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, fE1, fP1, fNP, energy);
  double I         = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, INorm, fNp2, energy);
  double ae        = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, fP2, fNp2, energy);
  double be        = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE3, fP3, fNp3, energy);
  pdf              = 0.0;
  if (I > 0.0) pdf = pe * sinh(sqrt(be * mid)) * exp(-mid / ae) / I;
  double pdfmid1   = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  if (fabs((pdf - pdfmid1) / pdfmid1) <= 5E-3) {
    return 0;
  }
  fEnergyFile5.push_back(mid);
  fEnergyPdfFile5.push_back(pdf);
  RecursionLinearFile5Watt(x1, mid, pdf1, pdf, energy);
  RecursionLinearFile5Watt(mid, x2, pdf, pdf2, energy);
  return 0;
}
