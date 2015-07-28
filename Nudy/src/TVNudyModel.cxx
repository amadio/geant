#include "TNudyCore.h"
#include "TVNudyModel.h"
// Header files required to display data
#include "TMath.h"
#include "TFile.h"
#include "TRandom3.h"
#include "Math/SpecFuncMathMore.h"
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TGraph2D.h>
#include <TStyle.h>
#include "base/Global.h"

using vecgeom::kPi;
using std::min;
using std::max;

ClassImp(TVNudyModel)

    //______________________________________________________________________________
    TVNudyModel::TVNudyModel() {
  fMAT = 0;
  fTemp = 0;
  fReaction = (Reaction_t)0;
  fMaterial = NULL;
  fProjectile = NULL;
  fE_file3 = NULL;
  fXSect_file3 = NULL;
  fEXSect_length = 0;
  f4nens = 0;
  fAPAlias = NULL;
  fEPtable = NULL;
  f5Tein = -1;
  f5Tel = 0;
  fPerc = NULL;
  fEPAlias = NULL;
  maxpop = 50;
  nens = 0;
  nperc = 25;
}

//_______________________________________________________________________________
TVNudyModel::TVNudyModel(TGeoElementRN *mat, Reaction_t reac, unsigned long temp, TParticlePDG *projectile,
                         TNudyEndfMat *material) {
  SetName(TNudyCore::Instance()->GetKey(mat, reac, temp));
  fEndf = mat->ENDFCode();
  fPdg = projectile->PdgCode();
  fMaterial = mat;
  fReaction = (Reaction_t)reac;
  fTemp = temp;
  fProjectile = projectile;
  // Make TGeoElementRN -> MAT map?
  fMAT = 0;
  fE_file3 = NULL;
  fXSect_file3 = NULL;
  fEXSect_length = 0;
  fEPtable = NULL;
  fPerc = NULL;
  fEPAlias = NULL;
  maxpop = 50;
  f5Tein = -1;
  f5Tel = 0;
  nens = 0;
  nperc = 25;
  if (material) {
    fMAT = material->GetMAT();
    ReadFile(material);
  }
}

//______________________________________________________________________________
TVNudyModel::~TVNudyModel() {
  if (fEPAlias) {
    delete[] fEPAlias;
    fEPAlias = 0;
  }
  if (fAPAlias) {
    delete[] fAPAlias;
    fAPAlias = 0;
  }
  if (fE_file3) {
    delete[] fE_file3;
    fE_file3 = 0;
  }
  if (fXSect_file3) {
    delete[] fXSect_file3;
    fXSect_file3 = 0;
  }
  if (fPerc) {
    delete[] fPerc;
    fPerc = 0;
  }
}

//______________________________________________________________________________
double TVNudyModel::GetEo(double ein) {
  int el, eh;
  if (f5Tein != ein) {
    Info("EIN", "%e --- %e", f5Tein, ein);
    f5Tein = ein;
    if (ein <= xengr[0]) {
      Warning("GetEo", "Incident energy %e below minimum threshold %e", ein, xengr[0]);
      el = 0;
    } else if (ein >= xengr[xengr.GetSize() - 1]) {
      Warning("GetEo", "Incident energy %e above maximum threshold %e", ein, xengr[xengr.GetSize() - 1]);
      el = xengr.GetSize() - 1;
    } else {
      el = TNudyCore::Instance()->BinarySearch(xengr.GetArray(), nens, ein);
      Info("Searching", "%d", el);
    }
    Info("", "%d - > %d", f5Tel, el);
    f5Tel = el;
  } else {
    el = f5Tel;
  }
  /*  eh = el+1;
    double pIndex = fRnd.Rndm();
    int pl = (pIndex * 25.0);
    int ph = pl + 1;
    return (TNudyCore::Instance()->BilinearInterploation(xengr[el],pl/25.0,xengr[eh],ph/25.0,
                                                         fPerc[el].GetAt(pl),fPerc[el].GetAt(ph),fPerc[eh].GetAt(pl),fPerc[eh].GetAt(ph),
                                                         ein,pIndex));
  */
  return fEPAlias[el].ImprovedInterpolation(f5Tein);
}
//_____________________________________________________________________________
double TVNudyModel::GetAo(double ein) {
  int el, eh;
  if (f4Tein != ein) {
    Info("EIN", "%e --- %e", f5Tein, ein);
    f4Tein = ein;
    if (ein <= fAPAlias[0].GetAlpha()) {
      Warning("GetAo", "Incident energy %e below minimum threshold %e", ein, fAPAlias[0].GetAlpha());
      el = 0;
    } else if (ein >= fAPAlias[f4nens - 1].GetAlpha()) {
      Warning("GetAo", "Incident energy %e above maximum threshold %e", ein, fAPAlias[f4nens - 1].GetAlpha());
      el = f4nens - 2;
    } else {
      el = TNudyCore::Instance()->BinarySearch(f4eins.GetArray(), f4nens, ein);
      Info("Searching", "%d", el);
    }
    f4Tel = el;
  } else {
    el = f4Tel;
  }
  return fAPAlias[el].ImprovedInterpolation(f4Tein);
}
//______________________________________________________________________________
TGeoElementRN *TVNudyModel::GetMaterial() {
  if (!fMaterial)
    fMaterial = TNudyCore::Instance()->GetMaterial(fEndf);
  return fMaterial;
}

//______________________________________________________________________________
const char *TVNudyModel::GetMaterialName() {
  if (!fMaterial)
    fMaterial = TNudyCore::Instance()->GetMaterial(fEndf);
  return fMaterial->GetName();
}

//______________________________________________________________________________
int TVNudyModel::GetZ() {
  if (!fMaterial)
    fMaterial = TNudyCore::Instance()->GetMaterial(fEndf);
  return fMaterial->Z();
}

//______________________________________________________________________________
int TVNudyModel::GetA() {
  if (!fMaterial)
    fMaterial = TNudyCore::Instance()->GetMaterial(fEndf);
  return (int)fMaterial->A();
}

//______________________________________________________________________________
int TVNudyModel::GetISO() {
  if (!fMaterial)
    fMaterial = TNudyCore::Instance()->GetMaterial(fEndf);
  return fMaterial->IsoNo();
}

//______________________________________________________________________________
Reaction_t TVNudyModel::GetReaction() {
  if (!fMaterial)
    fMaterial = TNudyCore::Instance()->GetMaterial(fEndf);
  return fReaction;
}

//______________________________________________________________________________
unsigned long TVNudyModel::GetTemp() {
  if (!fMaterial)
    fMaterial = TNudyCore::Instance()->GetMaterial(fEndf);
  return fTemp;
}

//______________________________________________________________________________
void TVNudyModel::ReadFile(TNudyEndfMat *material) {
  // Function to Read and Process File in ENDF Tape
  if (material->GetMAT() != fMAT) {
    printf("Wrong Material");
    return;
  }
  TIter iter(material->GetFiles());
  TNudyEndfFile *file;
  while ((file = (TNudyEndfFile *)iter.Next())) {
    // Read File data into class structures
    SetTitle(Form("%s:%d", GetTitle(), file->GetMF()));
    switch (file->GetMF()) {
    case 3:
      ReadFile3(file);
      break;
    case 4:
      ReadFile4(file);
      break;
    case 5:
      SetTitle(Form("%s(", GetTitle()));
      ReadFile5(file);
      SetTitle(Form("%s)", GetTitle()));
      break;
    }
  }
}

//______________________________________________________________________________
void TVNudyModel::ReadFile5(TNudyEndfFile *file) {
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;

  // nens = 0;
  // maxpop = 50;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    if (sec->GetMT() == (int)fReaction) {
      File5_Pass1(sec);
      File5_Pass2(sec);
      File5_Pass3();
    }
  }
}
//______________________________________________________________________________
void TVNudyModel::ReadFile4(TNudyEndfFile *file) {
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    if (sec->GetMT() == (int)fReaction) {
      TIter recIter(sec->GetRecords());
      TNudyEndfCont *header = (TNudyEndfCont *)recIter.Next();
      int LTT = sec->GetL2();
      int LI = header->GetL1();
      printf("LTT = %d LI = %d\n", LTT, LI);
      if (LTT == 1 && LI == 0) {
        TNudyEndfTab2 *subheader = (TNudyEndfTab2 *)recIter.Next();
        TArrayD ein(subheader->GetN2());
        TArrayD *lCoef = new TArrayD[subheader->GetN2()];
        for (int i = 0; i < subheader->GetN2(); i++) {
          TNudyEndfList *tab = (TNudyEndfList *)recIter.Next();
          lCoef[i].Set(tab->GetNPL());
          ein[i] = tab->GetC2();
          for (int j = 0; j < tab->GetNPL(); j++)
            lCoef[i].SetAt(tab->GetLIST(j), j);
        }
        for (int i = 0; i < ein.GetSize(); i++) {
          printf("Ein = %e\n", ein[i]);
          for (int j = 0; j < lCoef[i].GetSize(); j++) {
            printf("a%d = %e\n", j, lCoef[i].At(j));
          }
        }
      } else if (LTT == 3 && LI == 0) {
        TNudyEndfTab2 *lowE = (TNudyEndfTab2 *)recIter.Next();
        TArrayD ein(lowE->GetN2());
        TArrayD *lCoef = new TArrayD[lowE->GetN2()];
        for (int i = 0; i < lowE->GetN2(); i++) {
          TNudyEndfList *tab = (TNudyEndfList *)recIter.Next();
          lCoef[i].Set(tab->GetNPL());
          ein[i] = tab->GetC2();
          for (int j = 0; j < tab->GetNPL(); j++)
            lCoef[i].SetAt(tab->GetLIST(j), j);
        }
        TNudyEndfTab2 *highE = (TNudyEndfTab2 *)recIter.Next();
        f4nens = highE->GetN2();
        f4eins.Set(f4nens);
        fAPAlias = new TNudyAliasCont[f4nens];
        TArrayD *tVal = new TArrayD[highE->GetN2()];
        ein.Set(ein.GetSize() + highE->GetN2());
        for (int i = 0; i < highE->GetN2(); i++) {
          TNudyEndfTab1 *tab = (TNudyEndfTab1 *)recIter.Next();
          f4eins.SetAt(tab->GetC2(), i);
          ein[lowE->GetN2() + i] = tab->GetC2();
          tVal[i].Set(tab->GetNP() * 2);
          fAPAlias[i].Initialize(tab->Y(), tab->X(), tab->GetNP(), tab->GetC2(), i);
          for (int j = 0; j < tab->GetNP(); j++) {
            tVal[i].SetAt(tab->GetX(j), j * 2);
            tVal[i].SetAt(tab->GetY(j), j * 2 + 1);
          }
        }

        for (int i = 0; i < lowE->GetN2(); i++) {
          printf("Ein = %e\n", ein[i]);
          for (int j = 0; j < lCoef[i].GetSize(); j++) {
            printf("a%d = %e\n", j, lCoef[i].At(j));
          }
        }
        double h = 0;
        for (int i = 0; i < highE->GetN2(); i++) {
          printf("Ein = %e\n", ein[i + lowE->GetN2()]);
          h = 0;
          for (int j = 0; j < tVal[i].GetSize() / 2; j++) {
            if (j > 1)
              h += 0.5 * (tVal[i].At(j * 2) - tVal[i].At(j * 2 - 2)) * (tVal[i].At(j * 2 + 1) + tVal[i].At(j * 2 - 1));
            printf("X = %e Y = %e\n", tVal[i].At(j * 2), tVal[i].At(j * 2 + 1));
          }
          printf("Integral = %e", h);
        }
        TNudyAliasCont::BuildIntermediate(fAPAlias, f4nens);
      }
    }
  }
}
//______________________________________________________________________________
void TVNudyModel::File5_Pass3() {
  nperc = 25;
  fPerc = new TArrayD[nens];
  fEPAlias = new TNudyAliasCont[nens];
  for (int jen = 1; jen <= nens; jen++) {
    fPerc[jen - 1].Set(nperc + 1);
    int ibeg = 0, iend = 0;
    int ipoint = 0;
    double hintno = 0;
    double aX[nEout[jen - 1]];
    double aP[nEout[jen - 1]];
    for (int jp = 1; jp <= nEout[jen - 1] - 1; jp++) {
      hintno += 0.5 * (fEPtable[jen - 1].GetAt(2 * jp - 1) + fEPtable[jen - 1].GetAt(2 * jp + 1)) *
                (fEPtable[jen - 1].GetAt(2 * jp) - fEPtable[jen - 1].GetAt(2 * jp - 2));
    }
    if (fabs(hintno - 1) > 2e-2) {
      Warning("File5_Pass3", "Integral is not normalized at %e Integral = %e", xengr[jen - 1], hintno);
    }
    if (hintno <= 0) {
      for (int jp = 1; jp <= nperc + 1; jp++) {
        fPerc[jen - 1].SetAt(0, jp - 1);
      }
    } else {
      int np;
      for (np = 1; np <= nEout[jen - 1] - 1; np++) {
        if (fEPtable[jen - 1].GetAt(2 * np + 3) > 0) {
          //    Info("Pass3","First Value = %e E = %e",fEPtable[jen-1].GetAt(2*np+2),xengr[jen-1]);
          fPerc[jen - 1].SetAt(fEPtable[jen - 1].GetAt(2 * np), 0);
          ibeg = np;
          break;
        }
      }
      if (!(np <= nEout[jen - 1] - 1))
        Error("File5_Pass3", "No first point different from 0");
      for (np = nEout[jen - 1]; np >= 1; np--) {
        if (fEPtable[jen - 1].GetAt(2 * np - 1) > 0) {
          //  Info("Pass3","Last Value = %e E = %e",fEPtable[jen-1].GetAt(2*np-2),xengr[jen-1]);
          fPerc[jen - 1].SetAt(fEPtable[jen - 1].GetAt(2 * np - 2), nperc);
          iend = np;
          break;
        }
      }
      if (!(np >= 0))
        Error("File5_Pass3", "No last point different from 0");
      double hinteg = 0;
      double hintol = 0;
      ipoint = ibeg - 1;
      for (int jperc = 1; jperc <= nperc - 1; jperc++) {
        double percen = 1.0 * jperc / nperc;
        while (hinteg < percen && ipoint < iend) {
          hintol = hinteg;
          //      Info("F(x)","
          //  Info("int","ibeg = %d ipoint = %d, iend = %d , Int = %e",ibeg,ipoint,iend, hintol);
          ipoint++;
          hinteg += 0.5 * (fEPtable[jen - 1].GetAt(2 * ipoint - 1) + fEPtable[jen - 1].GetAt(2 * ipoint + 1)) *
                    (fEPtable[jen - 1].GetAt(2 * ipoint) - fEPtable[jen - 1].GetAt(2 * ipoint - 2)) / hintno;
        }
        double f1 = fEPtable[jen - 1].GetAt(2 * ipoint - 1) / hintno;
        double f2 = fEPtable[jen - 1].GetAt(2 * ipoint + 1) / hintno;
        double e1 = fEPtable[jen - 1].GetAt(2 * ipoint - 2);
        double e2 = fEPtable[jen - 1].GetAt(2 * ipoint);
        double delf = f2 - f1;
        double hintrp;
        if (delf != 0) {
          double dele = e2 - e1;
          double temp = 2 * (percen - hintol) * delf / dele;
          hintrp = e1 + dele * (sqrt(pow(f1, 2) + temp) - f1) / delf;
        } else {
          hintrp = e1 + (percen - hintol) / f1;
        }
        fPerc[jen - 1].SetAt(hintrp, jperc);
        //        printf("Eout - %e, P - %f\n",hintrp,percen);
      }
    }
    //    if(iend >= nEout[jen-1])
    //      iend = nEout[jen-1] - 1;
    //    while(iend > ibeg && fEPtable[jen-1].GetAt(2*iend+1) == 0) iend--;
    while (fEPtable[jen - 1].GetAt(2 * iend + 1) == 0)
      iend--;

    printf("ENERGY = %e\n", xengr[jen - 1]);
    for (int ja = 0; ja <= iend; ja++) {
      aX[ja] = fEPtable[jen - 1].GetAt(2 * ja);
      aP[ja] = fEPtable[jen - 1].GetAt(2 * ja + 1);
      printf("FX = %e FP = %e\n", aX[ja], aP[ja]);
    }
    fEPAlias[jen - 1].Initialize(aP, aX, iend + 1, xengr[jen - 1]);
  }
  TNudyAliasCont::BuildIntermediate(fEPAlias, nens);
}

//______________________________________________________________________________
void TVNudyModel::File5_Pass2(TNudyEndfSec *sec) {
  TIter recIter(sec->GetRecords());
  for (int k = 0; k < sec->GetN1(); k++) {
    TNudyEndfTab1 *header = (TNudyEndfTab1 *)recIter.Next();
    CheckLinear(header);
    double u = header->GetC1();
    SetTitle(Form("%s,%d", GetTitle(), header->GetL2()));
    if (header->GetL2() == 1) {
      TNudyEndfTab2 *record = (TNudyEndfTab2 *)recIter.Next();
      int nep = record->GetN2();
      double *engr = new double[nep];
      TList temp;
      for (int j = 1; j <= nep; j++) {
        TNudyEndfTab1 *row = (TNudyEndfTab1 *)recIter.Next();
        temp.Add(row);
        engr[j - 1] = row->GetC2();
      }
      for (int jg = 1; jg <= nens; jg++) {
        double energy = xengr[jg - 1];
        double ppe = TNudyCore::Instance()->Interpolate(header->NBT(), header->INT(), header->GetN1(), header->X(),
                                                        header->Y(), header->GetN2(), energy);
        if (ppe <= 0)
          continue;
        int ichan = TNudyCore::Instance()->BinarySearch(engr, nep, energy);
        if (engr[ichan] > energy || engr[ichan + 1] < energy) {
          Error("File5_Pass2", "Error in Interpolation %e does not lie between %e and %e", energy, engr[ichan],
                engr[ichan + 1]);
        }
        double xrat = (engr[ichan + 1] - energy) / (engr[ichan + 1] - engr[ichan]);
        double xra1 = 1 - xrat;
        if (min<double>(xrat, xra1) < 0 || max<double>(xrat, xra1) > 1) {
          Error("File5_Pass2", "Error in XRAT, XRA1 = %e, %e\nE = %e Ei = %e Eo = %e", xrat, xra1, energy, engr[ichan],
                engr[ichan + 1]);
        }
        double hint = 0;
        double eold = 0;
        double pold = 0;
        TNudyEndfTab1 *lfunc = (TNudyEndfTab1 *)temp.At(ichan);
        TNudyEndfTab1 *ufunc = (TNudyEndfTab1 *)temp.At(ichan + 1);
        for (int jef = 1; jef <= nEout[jg - 1]; jef++) {
          double efen = fEPtable[jg - 1].GetAt(2 * jef - 2);
          double pb1 = TNudyCore::Instance()->Interpolate(lfunc->NBT(), lfunc->INT(), lfunc->GetN1(), lfunc->X(),
                                                          lfunc->Y(), lfunc->GetN2(), efen);
          double pb2 = TNudyCore::Instance()->Interpolate(ufunc->NBT(), ufunc->INT(), ufunc->GetN1(), ufunc->X(),
                                                          ufunc->Y(), ufunc->GetN2(), efen);
          double prob = xrat * pb1 + xra1 * pb2;
          //    Info("File5_Pass2:Integration","e1,p1 = %e,%e , e2,p2 = %e,%e at E = %e Int = %e A=%e
          //    B=%e",eold,pold,efen,prob,energy, hint,pb1,pb2);
          fEPtable[jg - 1].GetArray()[2 * jef - 1] += ppe * prob;
          if (jef > 1) {
            hint = hint + (0.5 * (efen - eold) * (prob + pold));
          }
          eold = efen;
          pold = prob;
        }
        if (fabs(1 - hint) > 2e-2) {
          Error("File5_Pass2", "Integral error = %e at E = %e PPE = %e", hint * ppe, energy, ppe);
        } else {
          Info("File5_Pass2", "Succesful Integral %e at E = %e PPE = %e", hint, energy, ppe);
        }
      }
    } else if (header->GetL1() == 5) {
      double u = header->GetC1();
      TNudyEndfTab1 *temptab = (TNudyEndfTab1 *)recIter.Next();
      TNudyEndfTab1 *probtab = (TNudyEndfTab1 *)recIter.Next();
      for (int jg = 1; jg <= nens; jg++) {
        double energy = xengr[jg - 1];
        if (energy <= u)
          continue;
        double ppe = TNudyCore::Instance()->Interpolate(header->NBT(), header->INT(), header->GetN1(), header->X(),
                                                        header->Y(), header->GetN2(), energy);
        if (ppe <= 0)
          continue;
        double teta = TNudyCore::Instance()->Interpolate(temptab->NBT(), temptab->INT(), temptab->GetN1(), temptab->X(),
                                                         temptab->Y(), temptab->GetN2(), energy);
        double hint = 0;
        double pold = 0;
        double eold = 0;
        for (int jef = 1; jef <= nEout[jg - 1]; jef++) {
          double efen = fEPtable[jg - 1].GetAt(2 * jef - 2);
          if (!(efen >= (energy - u))) {
            double prob = TNudyCore::Instance()->Interpolate(probtab->NBT(), probtab->INT(), probtab->GetNR(),
                                                             probtab->X(), probtab->Y(), probtab->GetN2(), efen / teta);
            fEPtable[jg - 1].GetArray()[2 * jef - 2] += ppe * prob;
            if (jef > 1)
              hint += 0.5 * (efen - eold) * (prob + pold);
            eold = efen;
            pold = prob;
          } else {
            break;
          }
        }
        if (fabs(1 - hint) > 2e-2) {
          Error("File5_Pass2", "Integral error = %e at E = %e", hint * ppe, energy);
        } else {
          Info("File5_Pass2", "Succesful Integral %e at E = %e PPE = %e", hint, energy, ppe);
        }
      }
    } else if (header->GetL2() == 7) {
      double u = header->GetC1();
      TNudyEndfTab1 *temptab = (TNudyEndfTab1 *)recIter.Next();
      for (int jg = 1; jg <= nens; jg++) {
        double energy = xengr[jg - 1];
        double ppe = TNudyCore::Instance()->Interpolate(header->NBT(), header->INT(), header->GetN1(), header->X(),
                                                        header->Y(), header->GetN2(), energy);
        if (ppe <= 0)
          continue;
        double teta = TNudyCore::Instance()->Interpolate(temptab->NBT(), temptab->INT(), temptab->GetN1(), temptab->X(),
                                                         temptab->Y(), temptab->GetN2(), energy);
        double rede = (energy - u) / teta;
        double hnorm = pow(teta, 1.5) * (sqrt(kPi) / 2 * erf(sqrt(rede)) - sqrt(rede) * exp(-rede));
        double hint = 0;
        double pold = 0;
        double eold = 0;
        for (int jef = 1; jef <= nEout[jg - 1]; jef++) {
          double efen = fEPtable[jg - 1].GetAt(2 * jef - 2);
          if (!(efen >= (energy - u))) {
            double prob = sqrt(efen) * exp(-efen / teta) / hnorm;
            fEPtable[jg - 1].GetArray()[2 * jef - 1] += ppe * prob;
            if (jef > 1)
              hint += 0.5 * (efen - eold) * (prob + pold);
            eold = efen;
            pold = prob;
          } else {
            break;
          }
        }
        if (fabs(1 - hint) > 2e-2) {
          Error("File5_Pass2", "Integral error = %e at E = %e", hint * ppe, energy);
        } else {
          Info("File5_Pass2", "Succesful Integral %e at E = %e PPE = %e", hint, energy, ppe);
        }
      }
    } else if (header->GetL2() == 9) {
      double u = header->GetC1();
      TNudyEndfTab1 *temptab = (TNudyEndfTab1 *)recIter.Next();
      for (int jg = 1; jg <= nens; jg++) {
        double energy = xengr[jg - 1];
        if (energy <= u)
          continue;
        double ppe = TNudyCore::Instance()->Interpolate(header->NBT(), header->INT(), header->GetN1(), header->X(),
                                                        header->Y(), header->GetN2(), energy);
        if (ppe <= 0)
          continue;
        double teta = TNudyCore::Instance()->Interpolate(temptab->NBT(), temptab->INT(), temptab->GetN1(), temptab->X(),
                                                         temptab->Y(), temptab->GetN2(), energy);
        double rede = (energy - u) / teta;
        double hnorm;
        if (rede > 1e-6)
          hnorm = pow(teta, 2) * (1 - exp(-rede) * (1 + rede));
        else
          hnorm = pow(teta, 2) * 0.5 * pow(rede, 2) * (1 - rede);
        double hint = 0;
        double pold = 0;
        double eold = 0;
        for (int jef = 1; jef <= nEout[jg - 1]; jef++) {
          double efen = fEPtable[jg - 1].GetAt(2 * jef - 2);
          if (!(efen >= (energy - u))) {
            double prob = sqrt(efen) * exp(-efen / teta) / hnorm;
            fEPtable[jg - 1].GetArray()[2 * jef - 1] += ppe * prob;
            if (jef > 1)
              hint += 0.5 * (efen - eold) * (prob + pold);
            eold = efen;
            pold = prob;
          } else {
            break;
          }
        }
        if (fabs(1 - hint) > 2e-2) {
          Error("File5_Pass2", "Integral error = %e at E = %e", hint * ppe, energy);
        } else {
          Info("File5_Pass2", "Succesful Integral %e at E = %e PPE = %e", hint, energy, ppe);
        }
      }
    } else if (header->GetL2() == 11) {
      double u = header->GetC1();
      TNudyEndfTab1 *atab = (TNudyEndfTab1 *)recIter.Next();
      TNudyEndfTab1 *btab = (TNudyEndfTab1 *)recIter.Next();
      for (int jg = 1; jg <= nens; jg++) {
        double energy = xengr[jg - 1];
        if (energy <= u)
          continue;
        double ppe = TNudyCore::Instance()->Interpolate(header->NBT(), header->INT(), header->GetN1(), header->X(),
                                                        header->Y(), header->GetN2(), energy);
        if (ppe <= 0)
          continue;
        double a = TNudyCore::Instance()->Interpolate(atab->NBT(), atab->INT(), atab->GetN1(), atab->X(), atab->Y(),
                                                      atab->GetN2(), energy);
        double b = TNudyCore::Instance()->Interpolate(btab->NBT(), btab->INT(), btab->GetN1(), btab->X(), btab->Y(),
                                                      btab->GetN2(), energy);
        double ab4 = 0.25 * a * b;
        double sab4 = sqrt(ab4);
        double elim = energy - u;
        double rede = (elim) / a;
        double srede = sqrt(rede);
        double hnorm = 0.5 * a * sqrt(kPi) * sab4 * exp(ab4) * (erf(srede - sab4) + erf(srede + sab4)) -
                       a * exp(-rede) * sinh(sqrt(b * elim));
        double hint = 0;
        double pold = 0;
        double eold = 0;
        for (int jef = 1; jef <= nEout[jg - 1]; jef++) {
          double efen = fEPtable[jg - 1].GetAt(2 * jef - 2);
          if (!(efen >= (energy - u))) {
            double prob = exp(-efen / a) * sinh(sqrt(b * efen)) / hnorm;
            fEPtable[jg - 1].GetArray()[2 * jef - 1] += ppe * prob;
            if (jef > 1)
              hint += 0.5 * (efen - eold) * (prob + pold);
            eold = efen;
            pold = prob;
          } else {
            break;
          }
        }
        if (fabs(1 - hint) > 2e-2) {
          Error("File5_Pass2", "Integral error = %e at E = %e", hint * ppe, energy);
        } else {
          Info("File5_Pass2", "Succesful Integral %e at E = %e PPE = %e", hint, energy, ppe);
        }
      }
    } else if (header->GetL2() == 12) {
      double u = header->GetC1();
      TNudyEndfTab1 *temptab = (TNudyEndfTab1 *)recIter.Next();
      double ef[2] = {temptab->GetC1(), temptab->GetC2()};
      for (int jg = 1; jg <= nens; jg++) {
        double energy = xengr[jg - 1];
        if (energy <= u)
          continue;
        double elim = energy - u;
        double ppe = TNudyCore::Instance()->Interpolate(header->NBT(), header->INT(), header->GetN1(), header->X(),
                                                        header->Y(), header->GetN2(), energy);
        if (ppe <= 0)
          continue;
        double tm = TNudyCore::Instance()->Interpolate(temptab->NBT(), temptab->INT(), temptab->GetN1(), temptab->X(),
                                                       temptab->Y(), temptab->GetN2(), energy);
        double hint = 0;
        double *tmppb = new double[nEout[jg - 1]];
        for (int jef = 1; jef <= nEout[jg - 1]; jef++) {
          double efen = fEPtable[jg - 1].GetAt(2 * jef - 2);
          if (!(efen >= (elim))) {
            double prob = 0;
            for (int jhl = 0; jhl < 2; jhl++) {
              double u1 = pow(sqrt(efen) - sqrt(ef[jhl]), 2) / tm;
              double u2 = pow(sqrt(efen) + sqrt(ef[jhl]), 2) / tm;
              prob += 0.5 * ((pow(u2, 1.5) * ROOT::Math::expint(u2) + TMath::Gamma(1.5, u2)) -
                             (pow(u1, 1.5) * ROOT::Math::expint(u1) + TMath::Gamma(1.5, u1))) /
                      (3 * sqrt(ef[jhl] * tm));
            }
            tmppb[jef - 1] = prob;

            if (jef > 1)
              hint += 0.5 * (tmppb[jef - 2] + tmppb[jef - 1]) *
                      (fEPtable[jg - 1].GetAt(2 * jef - 1) - fEPtable[jg - 1].GetAt(2 * jef - 3));
          } else {
            break;
          }
        }
        for (int j = 1; j <= nEout[jg - 1]; j++)
          fEPtable[jg - 1].GetArray()[2 * j - 1] += ppe * tmppb[j - 1] / hint;
      }
    }
  }
}

//______________________________________________________________________________
void TVNudyModel::File5_Pass1(TNudyEndfSec *sec) {
  int index = 0;
  TIter recIter(sec->GetRecords());
  fEPtable = new TArrayD[500];
  for (int k = 0; k < sec->GetN1(); k++) {
    Info("error", "%d < %d", k, sec->GetN1());
    TNudyEndfTab1 *header = (TNudyEndfTab1 *)recIter.Next();
    CheckLinear(header);
    double u = header->GetC1();
    Info("File5_Pass1", "%s - File 5 - LF %d", GetName(), header->GetL2());
    if (header->GetL2() == 1) {
      TNudyEndfTab2 *range = (TNudyEndfTab2 *)recIter.Next();
      CheckLinear(range);
      int nz;
      TNudyEndfTab1 *row;
      for (nz = 0; nz < range->GetN2(); nz++) {
        row = (TNudyEndfTab1 *)recIter.Next();
        // Linearize Interplolation if it is histogram like
        Linearize(row);
        // Add upper threshold value
        /*	if(row->GetX(0) == 0 && row->GetY(0) == 0){
                  TNudyEndfTab1 *temp = new
           TNudyEndfTab1(row->GetC1(),row->GetC2(),row->GetL1(),row->GetL2(),row->GetN1(),row->GetN2()-1);
                  //memcpy(temp->X(),row->X()+1,sizeof(double)*row->GetN2()-1);
                  //		    memcpy(temp->Y(),row->Y()+1,sizeof(double)*row->GetN2()-1);
                  //		    memcpy(temp->INT(),row->INT(),sizeof(int)*row->GetN1());
                  for(int j = 1; j < row->GetN2(); j++){
                    temp->SetX(row->GetX(j),j-1);
                    temp->SetY(row->GetY(j),j-1);
                  }
                  for(int i = 0; i < row->GetN1(); i++){
                    temp->SetNBT(row->GetNBT(i)-1,i);
                    temp->SetINT(row->GetINT(i),i);
                  }
                  row->Equate(temp);
                  delete temp;
                }*/
        if (row->GetY(row->GetN2() - 1) > 0) {
          TNudyEndfTab1 *temp = new TNudyEndfTab1(row, row->GetN1(), row->GetN2() + 1);
          temp->SetX(temp->GetX(row->GetN2() - 1) * (1 + 1e-4), row->GetN2());
          temp->SetY(1e-10, row->GetN2());
          row->Equate(temp);
          delete temp;
        }
        double tempen = row->GetC2();
        double ppe = TNudyCore::Instance()->Interpolate(header->NBT(), header->INT(), header->GetNR(), header->X(),
                                                        header->Y(), header->GetN2(), tempen);
        if (ppe <= 0)
          continue;
        int exists = 0;
        for (int jen = 0; jen < nens; jen++) {
          if (tempen == xengr[jen]) {
            index = jen;
            exists = 1;
            break;
          }
        }
        if (!exists) {
          nens++;
          index = nens - 1;
          xengr.Set(nens);
          nEout.Set(nens);
          xengr[index] = tempen;
          Info("Pass1", "Adding %e energy at %d", tempen, index);
          nEout[index] = 0;
          double emin = max<double>(row->GetX(0) * (1 - 1e-5), 1e-5);
          double emax = row->GetX(row->GetN2() - 1) * (1 + 1e-5);
          double fact = exp(log(emax / emin) / (maxpop - 1));
          double ef = emin / fact;
          for (int jef = 0; jef < maxpop; jef++) {
            ef = ef * fact;
            if (!EoExists(index, ef)) {
              //		    Info("Pass1","Adding outgoing energy %e at %d",ef,index);
              nEout[index]++;
              fEPtable[index].Set(2 * nEout[index]);
              fEPtable[index].SetAt(ef, 2 * nEout[index] - 2);
              fEPtable[index].SetAt(0, 2 * nEout[index] - 1);
            }
          }
          for (int jeps = 1; jeps <= row->GetNP(); jeps++) {
            ef = row->GetX(jeps - 1);
            // Separating equal probability secondary energies
            if (jeps - 1 < row->GetNP() && ef == row->GetX(jeps)) {
              ef = ef * (1 - 5e-6);
              row->SetX(row->GetX(jeps) * (1 + 5e-6), jeps);
            }
            if (!EoExists(index, ef)) {
              nEout[index]++;
              fEPtable[index].Set(2 * nEout[index]);
              fEPtable[index].SetAt(ef, 2 * nEout[index] - 2);
              fEPtable[index].SetAt(0, 2 * nEout[index] - 1);
            }
          }
        }
      }
    } else if (header->GetL2() == 5) {
      TNudyEndfTab1 *theta = (TNudyEndfTab1 *)recIter.Next();
      FillGrid(u, theta->GetN2(), theta, header);
      TNudyEndfTab1 *dEdTheta = (TNudyEndfTab1 *)recIter.Next();
      for (index = 0; index < nens; index++) {
        double energy = xengr[index];
        double tet = TNudyCore::Instance()->Interpolate(theta->NBT(), theta->INT(), theta->GetN1(), theta->X(),
                                                        theta->Y(), theta->GetN2(), energy);
        for (int jef = 1; jef <= dEdTheta->GetN2(); jef++) {
          double ef = tet * dEdTheta->GetX(jef - 1);
          if (!EoExists(index, ef)) {
            nEout[index]++;
            fEPtable[index].Set(2 * nEout[index]);
            fEPtable[index].SetAt(ef, 2 * nEout[index] - 2);
            fEPtable[index].SetAt(0, 2 * nEout[index] - 1);
          }
        }
      }
    } else if (header->GetL2() == 7) {
      TNudyEndfTab1 *theta = (TNudyEndfTab1 *)recIter.Next();
      FillGrid(u, theta->GetN2(), theta, header);
    } else if (header->GetL2() == 9) {
      TNudyEndfTab1 *theta = (TNudyEndfTab1 *)recIter.Next();
      FillGrid(u, theta->GetN2(), theta, header);
    } else if (header->GetL2() == 11) {
      TNudyEndfTab1 *a = (TNudyEndfTab1 *)recIter.Next();
      FillGrid(u, a->GetN2(), a, header);
      TNudyEndfTab1 *b = (TNudyEndfTab1 *)recIter.Next();
      FillGrid(u, b->GetN2(), b, header);
    } else if (header->GetL2() == 12) {
      TNudyEndfTab1 *tm = (TNudyEndfTab1 *)recIter.Next();
      FillGrid(u, tm->GetN2(), tm, header);
    }
  }
  for (int i = 1; i <= nens - 1; i++) {
    for (int j = i + 1; j <= nens; j++) {
      if (xengr[i - 1] > xengr[j - 1]) {
        double temp = xengr[i - 1];
        xengr[i - 1] = xengr[j - 1];
        xengr[j - 1] = temp;
        int ntemp = nEout[i - 1];
        nEout[i - 1] = nEout[j - 1];
        nEout[j - 1] = ntemp;
        char buffer[sizeof(TArrayD)];
        memcpy(buffer, &fEPtable[i - 1], sizeof(TArrayD));
        memcpy(&fEPtable[i - 1], &fEPtable[j - 1], sizeof(TArrayD));
        memcpy(&fEPtable[j - 1], buffer, sizeof(TArrayD));
      }
    }
  }
  for (int jen = 0; jen < nens; jen++) {
    for (int j1 = 0; j1 < nEout[jen] - 1; j1++) {
      for (int j2 = j1 + 1; j2 < nEout[jen]; j2++) {
        if (fEPtable[jen].GetAt(2 * j1) > fEPtable[jen].GetAt(2 * j2)) {
          double tmp = fEPtable[jen].GetAt(2 * j1);
          fEPtable[jen].GetArray()[2 * j1] = fEPtable[jen].GetAt(2 * j2);
          fEPtable[jen].GetArray()[2 * j2] = tmp;
        }
      }
    }
  }
}

//______________________________________________________________________________
void TVNudyModel::PopulateGrid(int index) {
  int max = 50;
  double emin = 1e-5;
  double emax = 2e7;
  double fact, ef;
  if (emax > 0) {
    fact = exp(log(emax / emin) / (max - 1));
    ef = emin / fact;
  } else {
  }
  if (index > nEout.GetSize() - 1)
    nEout.Set(index + 1);
  nEout[index] = max;
  fEPtable[index].Set(2 * max);
  for (int jeps = 1; jeps <= max; jeps++) {
    ef = ef * fact;
    fEPtable[index].GetArray()[2 * jeps - 2] = ef;
    fEPtable[index].GetArray()[2 * jeps - 1] = 0;
  }
}

//______________________________________________________________________________
void TVNudyModel::FillGrid(double u, int nep, TNudyEndfTab1 *tab, TNudyEndfTab1 *pe) {
  int nene = 0;
  TArrayD tmpene;
  int j, keps, index;
  double ratmax = 2;
  double ratio, diff, ratiol, nadd, eadd, fact, ef, ethre1, ethre2;
  int maxene = 200;
  tmpene.Set(maxene);
  for (j = 1; j <= nep; j++) {
    double et = tab->GetX(j - 1);
    double ppe =
        TNudyCore::Instance()->Interpolate(pe->NBT(), pe->INT(), pe->GetNR(), pe->X(), pe->Y(), pe->GetN2(), et);
    if (ppe <= 0)
      continue;
    if (et - u < 1e-5) {
      et = max<double>(u * 1.001, u + 1e-4);
    }
    nene++;
    tmpene[nene - 1] = et;
  }
  int ntota = 0;
  for (j = 1; j <= nene - 1; j++) {
    ratio = tmpene[j] / tmpene[j - 1];
    diff = 2 * fabs(tab->GetY(j - 1) - tab->GetY(j)) / max<double>(1e-10, tab->GetY(j - 1) + tab->GetY(j));
    if (ratio > ratmax && diff > 1e-6) {
      ratiol = log(ratio);
      nadd = ratiol / log(ratmax) + 1;
      fact = exp(ratiol / nadd);
      eadd = tmpene[j - 1];
      for (int jad = 1; jad <= nadd - 1; jad++) {
        eadd = eadd * fact;
        tmpene[nene + ntota + jad - 2] = eadd;
      }
      ntota = ntota + nadd - 1;
    }
  }
  nene = nene + ntota;
  for (j = 1; j <= nene; j++) {
    int check = 0;
    for (int jen = 1; jen <= nens; jen++) {
      if (fabs(2 * (tmpene[j - 1] - xengr[jen - 1]) / (tmpene[j - 1] + xengr[jen - 1])) < 1e-7) {
        index = jen - 1;
        check = 1;
        break;
      }
    }
    if (!check) {
      nens++;
      index = nens - 1;
      xengr.Set(nens);
      xengr[index] = tmpene[j - 1];
      PopulateGrid(index);
    }
    double emin = 1e-5;
    double emax = max<double>(xengr[index] - u, emin + 1e-5);
    if (emax > 0) {
      fact = exp(log(emax / emin) / (maxene - 1));
      ef = emin / fact;
    } else {
      fact = 1;
      ef = 0;
    }
    for (int jeps = 1; jeps <= maxene; jeps++) {
      ef = ef * fact;
      if (!EoExists(index, ef)) {
        nEout[index]++;
        fEPtable[index].Set(2 * nEout[index]);
        fEPtable[index].SetAt(ef, 2 * nEout[index] - 2);
        fEPtable[index].SetAt(0, 2 * nEout[index] - 1);
      }
    }
    check = 1;
    ethre1 = emax * (1 - 1e-4);

    for (keps = 1; keps <= nEout[index]; keps++) {
      if (fEPtable[index].GetAt(2 * keps - 2) == ethre1) {
        check = 0;
        break;
      }
    }
    if (check) {
      nEout[index]++;
      fEPtable[index].Set(2 * nEout[index]);
      fEPtable[index].SetAt(ethre1, 2 * (nEout[index]) - 2);
      fEPtable[index].SetAt(0, 2 * (nEout[index]) - 1);
    }
    ethre2 = emax * (1 + 1e-4);
    check = 1;
    for (keps = 1; keps <= nEout[index]; keps++) {
      if (fEPtable[index].GetAt(2 * keps - 2) == ethre2) {
        check = 0;
        break;
      }
    }
    if (check) {
      nEout[index]++;
      fEPtable[index].Set(2 * nEout[index]);
      fEPtable[index].SetAt(ethre2, 2 * nEout[index] - 2);
      fEPtable[index].SetAt(0, 2 * nEout[index] - 1);
    }
  }
}

//______________________________________________________________________________
int TVNudyModel::EoExists(int index, double ef) {
  int check = 0;
  for (int keps = 1; keps <= nEout[index]; keps++) {
    if (fEPtable[index].GetAt(2 * keps - 2) == ef) {
      check = 1;
      break;
    }
  }
  return check;
}

//______________________________________________________________________________
void TVNudyModel::Linearize(TNudyEndfTab1 *tab) {
  if (!CheckLinear(tab)) {
    int islin = 0;
    int start = 1;

    int n2 = tab->GetN2();
    double *epval = new double[tab->GetN2() * 4];
    for (int i = 0; i < tab->GetN2(); i++) {
      epval[2 * i] = tab->GetX(i);
      epval[2 * i + 1] = tab->GetY(i);
    }
    for (int jr = 1; jr <= tab->GetN1(); jr++) {
      if (tab->GetINT(jr - 1) == 1) {
        islin = 1;
        for (int jp = start; jp <= tab->GetNBT(jr - 1) - 1; jp++) {
          n2 = n2 + 1;
          epval[2 * n2 - 2] = max<double>((1 - 1e-6) * epval[2 * jp], epval[jp * 2 - 2]);
          epval[2 * n2 - 1] = epval[jp * 2 - 1];
        }
      } else {
        start = tab->GetNBT(jr - 1);
        if (tab->GetINT(jr) != 2) {
          Warning("Linearize", "Unexpected interpolation law");
        }
      }
    }
    if (!islin)
      return;
    for (int j1 = 1; j1 <= n2; j1++) {
      for (int j2 = j1 + 1; j2 <= n2; j2++) {
        if (epval[2 * j1 - 2] > epval[2 * j2 - 2]) {
          double tmp1 = epval[2 * j1 - 2];
          double tmp2 = epval[2 * j1 - 1];
          epval[2 * j1 - 2] = epval[2 * j2 - 2];
          epval[2 * j1 - 1] = epval[2 * j2 - 1];
          epval[2 * j2 - 2] = tmp1;
          epval[2 * j2 - 1] = tmp2;
        }
      }
    }
    TNudyEndfTab1 *linear = new TNudyEndfTab1(tab->GetC1(), tab->GetC2(), tab->GetL1(), tab->GetL2(), 1, n2);
    linear->SetINT(2, 0);
    linear->SetNBT(n2, 0);
    for (int i = 0; i < n2; i++) {
      linear->SetX(epval[2 * i], i);
      linear->SetY(epval[2 * i + 1], i);
    }
    tab->Equate(linear);
    delete linear;
  }
}

//_______________________________________________________________________________
int TVNudyModel::CheckLinear(TNudyEndfTab1 *tab) {
  for (int nr = 0; nr < tab->GetNR(); nr++) {
    if (tab->GetINT(nr) != 2) {
      Error("CheckLinear", "Tab1 record data is not interpolated linearly");
      return 0;
    }
  }
  return 1;
}
int TVNudyModel::CheckLinear(TNudyEndfTab2 *tab) {
  for (int nr = 0; nr < tab->GetNR(); nr++) {
    if (tab->GetINT(nr) != 2) {
      Error("CheckLinear", "Tab2 record data is not interpolated linearly");
      return 0;
    }
  }
  return 1;
}

//_______________________________________________________________________________
void TVNudyModel::DisplayData(FileData_t file) {
  // This function displays data for a particular file stored in this model
  // Points drawn in graph are those from data.
  TCanvas *c1 = new TCanvas("c1", "Display Data", 200, 10, 700, 500);
  c1->SetFillColor(42);
  c1->SetGrid();
  TCanvas *canvas;
  TGraph2D *gr2;
  TGraph *gr = NULL;
  // Draw data
  switch (file) {
  case kReac_XSect:
    if (!fE_file3 || !fXSect_file3) {
      Error("DisplayData()", "This model does not contain any data about file 3");
      return;
    }
    c1->SetLogy();
    c1->SetLogx();
    gr = new TGraph(fEXSect_length, fE_file3, fXSect_file3);
    gr->SetTitle("File 3");
    gr->GetXaxis()->SetTitle("Energy");
    gr->GetYaxis()->SetTitle("Cross Section");
    break;
  case kEnergy_Dist:
    if (!fPerc) {
      Error("DisplayData()", "This model does not contain any data about file 5");
      return;
    }
    gr2 = new TGraph2D();
    for (int i = 0; i < nens; i++) {
      for (int j = 1; j <= 24; j++) {
        gr2->SetPoint(i * 24 + j, xengr[i], j / 24.0, fPerc[i].GetAt(j));
      }
    }
    gStyle->SetPalette(1);
    gr2->SetName(GetMaterialName());
    gr2->SetTitle("Energy Probability Distribution");
    gr2->GetXaxis()->SetTitle("Primary Energy (Ein)");
    gr2->GetYaxis()->SetTitle("Probability (P)");
    gr2->GetZaxis()->SetTitle("Secondary Energy (Eo)");
    gr2->Draw("surf1");
    break;
  default:
    break;
  }
  if (gr) {
    gr->SetLineColor(2);
    gr->SetLineWidth(1);
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    gr->Draw("al");
  }
  if (gr2) {
    gr2->Draw("surf1");
  }
  c1->Update();
  c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderSize(5);
  c1->Modified();
  c1->SaveAs(Form("%d_%d_%d_%d_%d.eps", fMAT, fTemp, fPdg, fReaction, file));
  c1->Close();
  if (!gr) {
    delete gr;
  }
  if (!gr2) {
    delete gr2;
  }
  delete c1;
}

//_______________________________________________________________________________
void TVNudyModel::DumpData(FileData_t file) {

  switch (file) {
  case kEnergy_Dist:
    printf("Energy Probability Distribution for %s\n\n", fMaterial->GetName());
    break;
  default:
    return;
  }
}

//______________________________________________________________________________
TArrayD *TVNudyModel::GetFile5Data(double ein) {
  int lo, hi, mid;
  int found = -1;
  lo = 0;
  hi = nens - 1;
  while (hi - lo > 1) {
    mid = (lo + hi) / 2;
    if (ein > xengr.GetAt(mid))
      lo = mid;
    else if (ein < xengr.GetAt(mid))
      hi = mid;
    else {
      found = mid;
      break;
    }
  }
  if (found >= 0)
    return &fEPtable[found];
  else {
    printf("hi=%d,lo=%d\n", hi, lo);
    Error("getFile5Data", "No such Ein in data");
    return NULL;
  }
}

//______________________________________________________________________________
TArrayD *TVNudyModel::GetFile5ProcessedData(double ein) {
  int lo, hi, mid;
  int found = -1;
  lo = 0;
  hi = nens - 1;
  while (hi - lo > 1) {
    mid = (lo + hi) / 2;
    if (ein > xengr.GetAt(mid))
      lo = mid;
    else if (ein < xengr.GetAt(mid))
      hi = mid;
    else {
      found = mid;
      break;
    }
  }
  if (found >= 0)
    return &fPerc[found];
  else {
    printf("hi=%d,lo=%d\n", hi, lo);
    Error("getFile5Data", "No such Ein in data");
    return NULL;
  }
}
