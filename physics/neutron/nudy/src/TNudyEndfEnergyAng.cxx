// 	This class is reconstructing probability tables for Energy distribution
// 	of the secondatries
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: March 22, 2016

#include "TList.h"
#include "Geant/TNudyEndfFile.h"
#include "Geant/TNudyEndfSec.h"
#include "Geant/TNudyEndfTab1.h"
#include "Geant/TNudyEndfTab2.h"
#include "Geant/TNudyEndfList.h"
#include "Geant/TNudyCore.h"
#include "Geant/TNudyEndfEnergyAng.h"
#include "Math/SpecFuncMathMore.h"

using namespace Nudy;
using namespace NudyPhysics;

#ifdef USE_ROOT
ClassImp(TNudyEndfEnergyAng)
#endif

TNudyEndfEnergyAng::TNudyEndfEnergyAng() : TNudyEndfSigma()
{
}

//______________________________________________________________________________
TNudyEndfEnergyAng::TNudyEndfEnergyAng(TNudyEndfFile *file, double iQValue[])
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    int NK = sec->GetN1();
    TIter recIter(sec->GetRecords());
    for (int k = 0; k < NK; k++) {
      TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
      div_t divr;
      double c[2];
      int nl[4];
      mtf[0]     = sec->GetMAT();
      mtf[1]     = sec->GetMT();
      mtf[2]     = sec->GetMF();
      AWRI       = sec->GetC2();
      divr       = div(sec->GetC1(), 1000);
      double ZA  = divr.quot;
      double AA  = divr.rem;
      double NA1 = AA - ZA;
      double AWP = tab1->GetC2();
      int ZAP    = tab1->GetC1();
      // int LIP = tab1->GetL1();
      int LAW = tab1->GetL2();
      // There is no structure for these laws
      if (LAW == 3 || LAW == 4 || LAW == 0) continue;
      fNR           = tab1->GetN1();
      fNP           = tab1->GetN2();
      divr          = div(ZAP, 1000);
      int particleA = divr.rem;
      int particleZ = divr.quot;
      int ID = 1;
      if (particleA == 1 && particleZ == 0) ID = 0;
      if (particleA == 4 && particleZ == 2) ID = 5;
      //---------------------------------------------
      switch (LAW) {
        case 1: 
        { 
          int NA, NEP, LANG;
          c[0]  = tab1->GetC1();
          c[1]  = tab1->GetC2();
          nl[0] = tab1->GetL1();
          nl[1] = 7;
          nl[2] = tab1->GetN1();
          nl[3] = tab1->GetN2();          
          // converting LAW 6 to LAW 7 after generating the distributions in LAB by multi-body kinematics 
          // given in LAW 6
          for (int i = 0; i < fNR; i++) {
            fNbt1.push_back(tab1->GetNBT(i));
            fInt1.push_back(tab1->GetINT(i));
          }
          for (int crs = 0; crs < fNP; crs++) {
            fE1.push_back(tab1->GetX(crs));
            fP1.push_back(tab1->GetY(crs));
          }
          tab1->SetCont(c[0],c[1],nl[0],nl[1],nl[2],nl[3]);
          for (int i = 0; i < fNR; i++) {
            tab1->SetNBT(fNbt1[i],i);
            tab1->SetINT(fInt1[i],i);
          }
          for (int crs = 0; crs < fNP; crs++) {
            tab1->SetX(fE1[crs],crs);
            tab1->SetY(fP1[crs],crs);
          }
          TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
          LANG                = tab2->GetL1();
          fNr2 = tab2->GetN1();
          fNp2 = tab2->GetN2();
          for (int lis = 0; lis < fNp2; lis++) {
            TNudyEndfList *header = (TNudyEndfList *)recIter.Next();
            fEin.push_back(header->GetC2());
            double sumein = 0;
            NA  = header->GetL2();
            NEP = header->GetN2();
            switch (LANG) {
              case 2:
              {
                switch (NA) {
                  case 1: // This law is tested with n-005_B_011.endf ENDF-B-VIII.0
                  {
                    for (int lis1 = 0; lis1 < NEP; lis1++) {
                      fEdes6.push_back(header->GetLIST(lis1 * 3 + 0));
                      fF06.push_back(header->GetLIST(lis1 * 3 + 1));
                      fR6.push_back(header->GetLIST(lis1 * 3 + 2));
                    }
                  }break;
                  case 2:
                  {
                    for (int lis1 = 0; lis1 < NEP; lis1++) {
                      fEdes6.push_back(header->GetLIST(lis1 * 4 + 0));
                      fF06.push_back(header->GetLIST(lis1 * 4 + 1));
                      fR6.push_back(header->GetLIST(lis1 * 4 + 2));
                      fA6.push_back(header->GetLIST(lis1 * 4 + 3));
                    }
                  }break;
                }
                double epsa = header->GetC2() * AWRI / (1. + AWRI);
                double AC = 1. + AA, ZC = ZA, NC = AC - ZC;
                double AB = AC - particleA, ZB = ZA - particleZ, NB = AB - ZB;
                double Ia = 0.0;
                double Ib = 0.0;
                double Sa = 15.68 * (AC - AA) - 28.07 * ((NC - ZC) * (NC - ZC) / AC - (NA1 - ZA) * (NA1 - ZA) / AA) -
                            18.56 * (pow(AC, 2/3) - pow(AA, 2/3)) +
                            33.22 * ((NC - ZC) * (NC - ZC) / pow(AC, 4/3) - (NA1 - ZA) * (NA1 - ZA) / pow(AA, 4/3)) -
                            0.717 * (ZC * ZC / pow(AC, 1/3) - ZA * ZA / pow(AA, 1/3)) +
                            1.211 * (ZC * ZC / AC - ZA * ZA / AA) - Ia;
                double Sb = 15.68 * (AC - AB) - 28.07 * ((NC - ZC) * (NC - ZC) / AC - (NB - ZB) * (NB - ZB) / AB) -
                            18.56 * (pow(AC, 2/3) - pow(AB, 2/3)) +
                            33.22 * ((NC - ZC) * (NC - ZC) / pow(AC, 4/3) - (NB - ZB) * (NB - ZB) / pow(AB, 4/3)) -
                            0.717 * (ZC * ZC / pow(AC, 1/3) - ZB * ZB / pow(AB, 1/3)) +
                            1.211 * (ZC * ZC / AC - ZB * ZB / AB) - Ib;
                double C1 = 0.04, C2 = 1.8E-6, C3 = 6.7E-7, Et1 = 130.0, Et3 = 41.0;
                double Mn1[6] = {1,1,1,1,1,0};
                double mn[6] = {0.5,1,1,1,1,2};
                int k1 = 0;
                do {
                  double x = -1. + k1 * 0.02;
                  double sumprob = 0;
                  for (int j = 0, Edes6Size = fEdes6.size(); j != Edes6Size; ++j) {
                    double epsb = fEdes6[j] * (AWP + AWRI) / (AWRI + 1. - AWP);
                    double ea   = epsa + Sa;
                    double eb   = epsb + Sb;
                    double R1   = std::min(ea, Et1), R3 = std::min(ea, Et3);
                    double X1   = R1 * eb / ea;
                    double X3   = R3 * eb / ea;
                    double a0   = C1 * X1 + C2 * X1 * X1 * X1 + C3 * Mn1[ID] * mn[ID] * pow(X3, 4);
                    double prob = (0.5 * a0 * fF06[j] / sinh(a0)) * (cosh(a0 * x) + fR6[j] * sinh(a0 * x));
                    if (NA == 2) prob = 
                      (0.5 * fA6[j] * fF06[j] / sinh(fA6[j])) * (cosh(fA6[j] * x) + fR6[j] * sinh(fA6[j] * x));
                    fEnergyFile5.push_back(fEdes6[j]);
                    fEnergyPdfFile5.push_back(prob);
                    sumprob += prob;
                  }
                  double fsum = 0.0;
                  for (int cr = 1, EnergyFile5Size = fEnergyFile5.size(); cr != EnergyFile5Size; ++cr) {
                    fsum +=
                        0.5 * (fEnergyPdfFile5[cr] + fEnergyPdfFile5[cr - 1]) * (fEnergyFile5[cr] - fEnergyFile5[cr - 1]);
                  }
                  sumein += fsum;
                  fCosIn.push_back(x);
                  fCosInPdf.push_back(fsum);
                  for (int i = 0, EnergyFile5Size = fEnergyFile5.size(); i != EnergyFile5Size; ++i) {
                    fEoute.push_back(fEnergyFile5[i]);
                    fPdfe.push_back(fEnergyPdfFile5[i]);
                  }
                  fEout2de.push_back(fEoute);
                  fPdf2de.push_back(fPdfe);
                  fEnergyFile5.clear();
                  fEnergyPdfFile5.clear();
                  fEoute.clear();
                  fPdfe.clear();
                  k1++;
                } while (k1 < 101);
                // creating records similar to law 7 
                TNudyEndfTab2 *secTab21 = new TNudyEndfTab2();
                int NMU = fCosIn.size();
                c[0]  = 0;
                c[1]  = fEin[lis];
                nl[0] = 0;
                nl[1] = 0;
                nl[2] = 1;
                nl[3] = NMU;          
                secTab21->SetCont(c[0],c[1],nl[0],1,nl[2],nl[3]);
                secTab21->SetNBT(NMU,0);
                secTab21->SetINT(2,0);
                sec->AddBefore(header,secTab21);
                
                TNudyEndfTab1 *secTab1(new TNudyEndfTab1[NMU]);
                for (int j = 0; j != NMU; ++j) {
                  CreateTab1(&secTab1[j], fEout2de[j], fPdf2de[j], mtf, 2, fCosIn[j]);
                  if (j == 0) {
                    sec->AddAfter(secTab21,&secTab1[j]);
                  } else {
                    sec->AddAfter(&secTab1[j-1],&secTab1[j]); 
                  }
                }
                sec->RemoveObj(header);
                fEout2de.clear();
                fPdf2de.clear();
                fCosIn.clear();
                fCosInPdf.clear();
                fEdes6.clear();
                fF06.clear();
                fR6.clear();
                fA6.clear();
              }break;
              case 1: // this law tested with n-005_B_010.endf file
              { 
                for (int lis1 = 0; lis1 < NEP; lis1++) {
                  fEdes6.push_back(header->GetLIST(lis1 * (NA + 2) + 0));
                  for (int j = 0; j < NA +1; j++) {
                    fLegendCoef.push_back(header->GetLIST(lis1 * (NA + 2) + j + 1));
                  }
                  fLegendCoef1.push_back(fLegendCoef);
                  fLegendCoef.clear();
                }
                int k1     = 0;
                double fme = 0;
                do {
                  double sumprob = 0;
                  double x       = -1. + k1 * 0.02;
                  for (int m = 0, Edes6Size = fEdes6.size(); m != Edes6Size; ++m) {
                    fme = 0.5;
                    for (int j = 0, LegendCoefSize = fLegendCoef1[m].size(); j != LegendCoefSize; ++j) {
                      double leg = ROOT::Math::legendre(j + 1, x);
                      fme       += 0.5 * (2. * (j + 1) + 1.) * fLegendCoef1[m][j] * leg;
                    }
                    if (fme < 0) fme = 0.5;
                    if (fme > 1E-15) {
                      fEnergyFile5.push_back(fEdes6[m]);
                      fEnergyPdfFile5.push_back(fme);
                      sumprob += fme;
                    }
                  }
                  double fsum = 0.0;
                  for (int cr = 1, EnergyFile5Size = fEnergyFile5.size(); cr != EnergyFile5Size; ++cr) {
                    fsum +=
                    0.5 * (fEnergyPdfFile5[cr] + fEnergyPdfFile5[cr - 1]) * (fEnergyFile5[cr] - fEnergyFile5[cr - 1]);
                  }
                  sumein += fsum;
                  fCosIn.push_back(x);
                  fCosInPdf.push_back(fsum);
                  for (int i = 0, EnergyFile5Size = fEnergyFile5.size(); i != EnergyFile5Size; ++i) {
                    fEoute.push_back(fEnergyFile5[i]);
                    fPdfe.push_back(fEnergyPdfFile5[i]);
                  }
                  fEout2de.push_back(fEoute);
                  fPdf2de.push_back(fPdfe);

                  fEnergyFile5.clear();
                  fEnergyPdfFile5.clear();
                  fEoute.clear();
                  fPdfe.clear();
                  k1++;
                } while (k1 < 101);
                // creating records similar to law 7 
                TNudyEndfTab2 *secTab21 = new TNudyEndfTab2();
                int NMU = fCosIn.size();
                c[0]  = 0;
                c[1]  = fEin[lis];
                nl[0] = 0;
                nl[1] = 0;
                nl[2] = 1;
                nl[3] = NMU;          
                secTab21->SetCont(c[0],c[1],nl[0],1,nl[2],nl[3]);
                secTab21->SetNBT(NMU,0);
                secTab21->SetINT(2,0);
                sec->AddBefore(header,secTab21);
                
                TNudyEndfTab1 *secTab1(new TNudyEndfTab1[NMU]);
                for (int j = 0; j != NMU; ++j) {
                  CreateTab1(&secTab1[j], fEout2de[j], fPdf2de[j], mtf, 2, fCosIn[j]);
                  if (j == 0) {
                    sec->AddAfter(secTab21,&secTab1[j]);
                  } else {
                    sec->AddAfter(&secTab1[j-1],&secTab1[j]); 
                  }
                }
                sec->RemoveObj(header);
                fEout2de.clear();
                fPdf2de.clear();
                fCosIn.clear();
                fCosInPdf.clear();
                fEdes6.clear();
                fLegendCoef1.clear();
              }break;
            } 
          }
          fEin.clear();
          fNbt1.clear();
          fInt1.clear();
          fE1.clear();
          fP1.clear();
        } break;
        case 2: // this law testedwith n-005_B_010.endf file
        { 
          ProcessTab1(tab1, fNR, fNP, fNbt1, fInt1, fEint, fYi);
          fNbt1.clear();
          fInt1.clear();          
          TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
          ProcessTab2(tab2, fNr2, fNp2, fNbt2, fInt2);
          fNbt2.clear();
          fInt2.clear();          
          for (int lis = 0; lis < fNp2; lis++) {
            TNudyEndfList *header = (TNudyEndfList *)recIter.Next();
            fEin.push_back(header->GetC2());
            int LANG = header->GetL1();
            fNP      = header->GetN2();
            switch (LANG) {
              case 0:
              {
                for (int j = 0; j < fNP; j++) {
                  fLegendCoef.push_back(header->GetLIST(j));
                }
                LegendreToCosine(fCosFile4, fCosPdfFile4);
                for(int l = 0; l < fNP -1; l++){
                  RecursionLinearLeg1D(fCosFile4[l], fCosFile4[l+1], fCosPdfFile4[l], fCosPdfFile4[l+1]);
                }
                TNudyCore::Instance()->Sort(fCosFile4, fCosPdfFile4);
                TNudyEndfTab1 *secTab1 = new TNudyEndfTab1();
                CreateTab1(secTab1, fCosFile4, fCosPdfFile4, mtf, 2, fEin[lis]);
                sec->AddBefore(header,secTab1);
                sec->RemoveObj(header);                
                fCosFile4.clear();
                fCosPdfFile4.clear();
                fLegendCoef.clear();
              }break;
              default:
              {
                for (int i = 0; i < fNP; i++) {
                  fCosFile4.push_back(header->GetLIST(2 * i + 0));
                  fCosPdfFile4.push_back(header->GetLIST(2 * i + 1));
                }
                int size = fCosFile4.size();
                fNR      = 1;
                fNbt1.push_back(size);
                LANG = LANG % 10;
                fInt1.push_back(LANG);
                for (int l = 0; l < size - 1; l++) {
                  if (LANG == 1) continue;
                    RecursionLinearProb(fCosFile4[l], fCosFile4[l + 1], fCosPdfFile4[l], fCosPdfFile4[l + 1]);
                }
                TNudyCore::Instance()->Sort(fCosFile4, fCosPdfFile4);
                TNudyEndfTab1 *secTab1 = new TNudyEndfTab1();
                CreateTab1(secTab1, fCosFile4, fCosPdfFile4, mtf, 2, fEin[lis]);
                sec->AddBefore(header,secTab1);
                sec->RemoveObj(header);                
                fCosFile4.clear();
                fCosPdfFile4.clear();
                fNbt1.clear();
                fInt1.clear();
              }break;
            }
          }
          fEin.clear();
          fNbt2.clear();
          fInt2.clear();
        } break;
        case 3:
        case 4: 
        {
          for (int cr = 0; cr < fNR; cr++) {
            fNbt1.push_back(tab1->GetNBT(cr));
            fInt1.push_back(tab1->GetINT(cr));
          }
          for (int crs = 0; crs < fNP; crs++) {
            fE1.push_back(tab1->GetX(crs));
            fP1.push_back(tab1->GetY(crs));
            fEin.push_back(tab1->GetX(crs));
          }
          fNbt1.clear();
          fInt1.clear();
          fE1.clear();
          fP1.clear();
          fEin.clear();
          //---------------------------------------------
        } break;
        case 5: // charge particle scattering law is not implimented yet
        {
        } break;
        case 6: // this law is tested for n-001_H_002.endf file
        {  
          int NPSX;
          double APSX;
          c[0]  = tab1->GetC1();
          c[1]  = tab1->GetC2();
          nl[0] = tab1->GetL1();
          nl[1] = 7;
          nl[2] = tab1->GetN1();
          nl[3] = tab1->GetN2();          
          // converting LAW 6 to LAW 7 after generating the distributions in LAB by multi-body kinematics 
          // given in LAW 6
          for (int i = 0; i < fNR; i++) {
            fNbt1.push_back(tab1->GetNBT(i));
            fInt1.push_back(tab1->GetINT(i));
          }
          for (int crs = 0; crs < fNP; crs++) {
            fE1.push_back(tab1->GetX(crs));
            fP1.push_back(tab1->GetY(crs));
          }
          tab1->SetCont(c[0],c[1],nl[0],nl[1],nl[2],nl[3]);
          for (int i = 0; i < fNR; i++) {
            tab1->SetNBT(fNbt1[i],i);
            tab1->SetINT(fInt1[i],i);
          }
          for (int crs = 0; crs < fNP; crs++) {
            tab1->SetX(fE1[crs],crs);
            tab1->SetY(fP1[crs],crs);
          }
          TNudyEndfCont *cont   = (TNudyEndfCont *)recIter.Next();
          APSX                  = cont->GetC1();
          NPSX                  = cont->GetN2();
          double energy         = std::fabs(-(1. + AWRI) * iQValue[sec->GetMT()] / AWRI);
          double eout           = 1E-5;
          energy               += 0.001 * (energy);
          do {
            double probEn = 0;
            double Ea     = AWRI * energy / (1. + AWRI) + iQValue[sec->GetMT()];
            double EiMax  = (APSX - AWP) * Ea / APSX;
            double EStar  = energy *AWP/((AWRI+1)*(AWRI+1));
            double C[3]   = {4. / (PI * EiMax * EiMax), 105. / (32. * pow(EiMax, 3.5)), 256. / (14. * PI * pow(EiMax, 5))};
            int k1        = 0;
            double rn1 = 0, rn2 = 0, rn3 = 0, rn4 = 0, rn5 = 0, rn6 = 0, rn7 = 0, rn8 = 0, rn9 = 0;
            double p = 0, xp = 0, yp = 0, tp = 0, rn12 = 0, rn34 = 0;
            // sampling of the outgoing energy is sampled based on R28 of LA-9721-MS.
            do {
              double cosCM = -1 + k1 * 0.02;
              int iei = 0;
              double probCos = 0;
              do {
                do {
                  rn1  = fRng.uniform();
                  rn2  = fRng.uniform();
                  rn12 = rn1 * rn1 + rn2 * rn2;
                } while (rn12 > 1);
                do {
                  rn3  = fRng.uniform();
                  rn4  = fRng.uniform();
                  rn34 = rn3 * rn3 + rn4 * rn4;
                } while (rn34 > 1);
                rn5              = fRng.uniform();
                rn6              = fRng.uniform();
                rn7              = fRng.uniform();
                rn8              = fRng.uniform();
                if (NPSX == 3) p = rn5;
                if (NPSX == 4) p = rn5 * rn6;
                if (NPSX == 5) p = rn5 * rn6 * rn7 * rn8;
                rn9              = fRng.uniform();
                xp               = -rn1 * log(rn12) / rn12 - log(rn9);
                yp               = -rn3 * log(rn34) / rn34 - log(p);
                tp               = xp / (xp + yp);
                eout             = tp * EiMax;
                // double ppeCm     = 0.0;
                // ppeCm            = 2 * C[NPSX - 3] * sqrt(eout) * pow((EiMax - eout), 1.5 * NPSX - 4);
                double ppeLab    = 0.0;
                double EMaxLimit = EiMax - (EStar + eout -2 * cosCM * sqrt(EStar * eout));
                if (EMaxLimit > 0) {
                  ppeLab           = 2 * C[NPSX - 3] * sqrt(eout) * pow(EMaxLimit, 1.5 * NPSX - 4);                
                  fEnergyFile5.push_back(eout);
                  fEnergyPdfFile5.push_back(ppeLab);
                  probCos += ppeLab;
                }
                iei++;
              } while (iei < 200);
              if (probCos > 0) {
                probEn += probCos;
                TNudyCore::Instance()->Sort(fEnergyFile5, fEnergyPdfFile5);
                fCosIn.push_back(cosCM);
                fCosInPdf.push_back(probCos);
                for (int i = 0, EnergyFile5Size = fEnergyFile5.size(); i != EnergyFile5Size; ++i) {
                  fEoute.push_back(fEnergyFile5[i]);
                  fPdfe.push_back(fEnergyPdfFile5[i]);
                }
                fEout2de.push_back(fEoute);
                fPdf2de.push_back(fPdfe);
              }
              fEnergyFile5.clear();
              fEnergyPdfFile5.clear();
              fEoute.clear();
              fPdfe.clear();
              k1++;
            } while (k1 < 101);
            if (probEn > 0) {
              fEin.push_back(energy);
              fCos2d.push_back(fCosIn);
              fCosinpdf2d.push_back(fCosInPdf);
              fEout3de.push_back(fEout2de);
              fPdf3de.push_back(fPdf2de);
            }
            fEout2de.clear();
            fPdf2de.clear();
            fCosIn.clear();
            fCosInPdf.clear();
            energy *= 2;
          } while (energy < fE1[fNP - 1] + iQValue[sec->GetMT()]);
          TNudyEndfTab2 *secTab2 = new TNudyEndfTab2();
          int NE = fEin.size();
          CreateTab2(secTab2,NE);
          sec->AddBefore(cont,secTab2);
          sec->RemoveObj(cont);
          TNudyEndfTab2 *secTab21(new TNudyEndfTab2[NE]);
          int imax = 0;
          TNudyEndfTab1 *secTab1(new TNudyEndfTab1[NE*101]);
          for (int i = 0; i != NE; ++i) {
            int NMU = fCos2d[i].size();
            c[0]  = 0;
            c[1]  = fEin[i];
            nl[0] = 0;
            nl[1] = 0;
            nl[2] = 1;
            nl[3] = NMU;          
            secTab21[i].SetCont(c[0],c[1],nl[0],1,nl[2],nl[3]);
            secTab21[i].SetNBT(NMU,0);
            secTab21[i].SetINT(2,0);
            if (i == 0) {
              sec->AddAfter(secTab2,&secTab21[i]);
            } else {
              sec->AddAfter(&secTab1[imax-1],&secTab21[i]); 
            }
            for (int j = 0; j != NMU; ++j) {
              CreateTab1(&secTab1[imax + j], fEout3de[i][j], fPdf3de[i][j], mtf, 2, fCos2d[i][j]);
              if (j == 0) {
                sec->AddAfter(&secTab21[i],&secTab1[imax + j]);
              } else {
                sec->AddAfter(&secTab1[imax + j - 1], &secTab1[imax + j]); 
              }
            }
            imax += NMU;
          }
          fEin.clear();
          fEout3de.clear();
          fPdf3de.clear();
          fCos2d.clear();
          fCosinpdf2d.clear();
          fNbt1.clear();
          fInt1.clear();
          fE1.clear();
          fP1.clear();
        } break;
        case 7: // laboratory angle and energy law
        {
          ProcessTab1(tab1, fNR, fNP, fNbt1, fInt1, fE1, fP1);
          TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
          ProcessTab2(tab2, fNr2, fNp2, fNbt2, fInt2);
          for (int cr1 = 0; cr1 < fNp2; cr1++) {
            TNudyEndfTab2 *tab3 = (TNudyEndfTab2 *)recIter.Next();
            fEin.push_back(tab3->GetC2());
            int NMU = tab3->GetN2();
            for (int cr2 = 0; cr2 < NMU; cr2++) {
              TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
              ProcessTab1(tab12, fNr3, fNp3, fNbt3, fInt3, fEnergyFile5, fEnergyPdfFile5);
              double fsum = 0.0;
              for (int cr = 1; cr < fNp3; cr++) {
                fsum += 0.5 * (fEnergyPdfFile5[cr] + fEnergyPdfFile5[cr - 1]) * (fEnergyFile5[cr] - fEnergyFile5[cr - 1]);
              }
              fCosInPdf.push_back(fsum);
              int law = 2;
              for (int cr = 0; cr < fNp3 - 1; cr++) {
                for (int i = 0; i < fNr3; i++) {
                  if (fNbt3[i] > cr) {
                    law = fInt3[i];
                    break;
                  }
                }
                if (law == 1) continue;                
                RecursionLinearFile5Prob(fEnergyFile5[cr], fEnergyFile5[cr + 1], fEnergyPdfFile5[cr],
                                        fEnergyPdfFile5[cr + 1]);
              }
              TNudyCore::Instance()->Sort(fEnergyFile5, fEnergyPdfFile5);
              ModifyTab1(tab12, fEnergyFile5, fEnergyPdfFile5, fEin[cr1]);
              fEnergyFile5.clear();
              fEnergyPdfFile5.clear();
              fNbt3.clear();
              fInt3.clear();
            }
          }
          fEin.clear();
          fNbt1.clear();
          fInt1.clear();
          fE1.clear();
          fP1.clear();
          fNbt2.clear();
          fInt2.clear();
        } break;
      }
    }
  }
}
//------------------------------------------------------------------------------------------------------
TNudyEndfEnergyAng::~TNudyEndfEnergyAng()
{
  fCosFile4.shrink_to_fit();
  fEnergyFile5.shrink_to_fit();
  fCosPdfFile4.shrink_to_fit();
  fEnergyPdfFile5.shrink_to_fit();
  fEnergyCdfFile5.shrink_to_fit();
  fEdes6.shrink_to_fit();
  fF06.shrink_to_fit();
  fR6.shrink_to_fit();
  fA6.shrink_to_fit();
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
  fEin.shrink_to_fit();
  fLegendCoef1.shrink_to_fit();
  fCosIn.shrink_to_fit();
  fCosInPdf.shrink_to_fit();
  fEoute.shrink_to_fit();
  fPdfe.shrink_to_fit();
  fCos2d.shrink_to_fit();
  fLegendCoef.shrink_to_fit();
  fEout2de.shrink_to_fit();
  fPdf2de.shrink_to_fit();
  fEout3de.shrink_to_fit();
  fPdf3de.shrink_to_fit();
}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfEnergyAng::ProcessTab2(Nudy::TNudyEndfTab2 *tab2, int &NR,int &NP,rowint &fnbt, rowint &fint)
{
  NudyPhysics::TNudyEndfSigma::ProcessTab2(tab2, NR, NP, fnbt, fint);
}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfEnergyAng::ProcessTab1(Nudy::TNudyEndfTab1 *tab1,int &NR,int &NP,rowint &fnbt, 
                                  rowint &fint,rowd &x1, rowd &x2)
{  
  NudyPhysics::TNudyEndfSigma::ProcessTab1(tab1, NR, NP, fnbt, fint, x1, x2);
}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfEnergyAng::ModifyTab1(TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, double &x)
{
  NudyPhysics::TNudyEndfSigma::ModifyTab1(secTab1, x1, x2, x);
}  
// -------------------------------------------------------------------------------------------------------
void TNudyEndfEnergyAng::CreateTab2(TNudyEndfTab2 *secTab2, int &NE)
{
  NudyPhysics::TNudyEndfSigma::CreateTab2(secTab2, NE);
}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfEnergyAng::CreateTab1(TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, int mtf[], int law, double &x)
{
  NudyPhysics::TNudyEndfSigma::CreateTab1(secTab1, x1, x2, mtf, law, x);
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfEnergyAng::RecursionLinearProb(double x1, double x2, double pdf1, double pdf2)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2) return 0;
  pdf            = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, fCosFile4, fCosPdfFile4, fNP, mid);
  double pdfmid1 = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  if (fabs((pdf - pdfmid1) / pdfmid1) <= 1E-3) {
    return 0;
  }
  fCosFile4.push_back(mid);
  fCosPdfFile4.push_back(pdf);
  RecursionLinearProb(x1, mid, pdf1, pdf);
  RecursionLinearProb(mid, x2, pdf, pdf2);
  return 0;
}

//------------------------------------------------------------------------------------------------------
double TNudyEndfEnergyAng::RecursionLinearFile5Prob(double x1, double x2, double pdf1, double pdf2)
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
//------------------------------------------------------------------------------------------------------
double TNudyEndfEnergyAng::RecursionLinearLeg1D(double x1, double x2, double pdf1, double pdf2)
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
//------------------------------------------------------------------------------------------------------
void TNudyEndfEnergyAng::LegendreToCosine(rowd &x1, rowd &x2)
{
  int k1     = 0;
  double fme = 0.0;
  do {
    fme      = 1.0;
    double x = -1. + k1 * 0.02;
    for (unsigned long j = 0; j < fLegendCoef.size(); j++) {
      double leg = ROOT::Math::legendre(j + 1, x);
      fme       += 0.5 * (2. * (j + 1) + 1.) * fLegendCoef[j] * leg;
    }
    if (fme > 0.0) {
      x1.push_back(x);
      x2.push_back(fme);
    }
    k1++;
  } while (k1 < 101);
}

