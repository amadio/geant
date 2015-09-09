#include "TPFstate.h"
#include "TFinState.h"
#include <TFile.h>
#include <TMath.h>
#include <TRandom.h>

using std::max;

int TPFstate::fVerbose = 0;

ClassImp(TPFstate)

    //_________________________________________________________________________
    TPFstate::TPFstate()
    : fPDG(0), fNEbins(0), fNReac(0), fNEFstat(0), fNFstat(0), fEmin(0), fEmax(0), fEilDelta(0),
      fEGrid(TPartIndex::I()->EGrid()), fFstat(0), fRestCaptFstat(0) {
  int np = TPartIndex::I()->NProc();
  while (np--)
    fRdict[np] = fRmap[np] = -1;
}

//_________________________________________________________________________
TPFstate::TPFstate(int pdg, int nfstat, int nreac, const int dict[])
    : fPDG(pdg), fNEbins(TPartIndex::I()->NEbins()), fNReac(nreac), fNEFstat(nfstat), fNFstat(fNEbins * fNEFstat),
      fEmin(TPartIndex::I()->Emin()), fEmax(TPartIndex::I()->Emax()), fEilDelta((fNEbins - 1) / log(fEmax / fEmin)),
      fEGrid(TPartIndex::I()->EGrid()), fFstat(new TFinState[fNFstat]), fRestCaptFstat(0) {
  int np = TPartIndex::I()->NProc();
  while (np--) {
    fRdict[dict[np]] = np;
    fRmap[np] = dict[np];
  }
  // consistency
  for (int i = 0; i < fNReac; ++i)
    if (fRdict[fRmap[i]] != i)
      Fatal("SetPartXS", "Dictionary mismatch for!");
}

//_________________________________________________________________________
TPFstate::~TPFstate() {
  delete[] fFstat;
  if (fRestCaptFstat)
    delete fRestCaptFstat;
}

//_________________________________________________________________________
void TPFstate::SetRestCaptFstate(const TFinState &finstate) {
  if (finstate.GetNsecs()) {
    fRestCaptFstat = new TFinState();
    fRestCaptFstat->SetFinState(finstate);
  } else {
    fRestCaptFstat = 0;
  }
}

//_________________________________________________________________________
bool TPFstate::SetPart(int pdg, int nfstat, int nreac, const int dict[]) {
  fPDG = pdg;
  fNEbins = TPartIndex::I()->NEbins();
  fNReac = nreac;
  fNEFstat = nfstat;
  fNFstat = fNEbins * fNReac;
  fEmin = TPartIndex::I()->Emin();
  fEmax = TPartIndex::I()->Emax();
  fEGrid = TPartIndex::I()->EGrid();
  fEilDelta = (fNEbins - 1) / log(fEmax / fEmin);
  fFstat = new TFinState[fNFstat];

  int np = TPartIndex::I()->NProc();
  while (np--) {
    fRdict[dict[np]] = np;
    fRmap[np] = dict[np];
  }
  // consistency
  for (int i = 0; i < fNReac; ++i)
    if (fRdict[fRmap[i]] != i)
      Fatal("SetPart", "Dictionary mismatch for!");
  return kTRUE;
}

//_________________________________________________________________________
bool TPFstate::SetPart(int pdg, int nfstat, int nreac, const int dict[], TFinState vecfs[]) {
  fPDG = pdg;
  fNEbins = TPartIndex::I()->NEbins();
  fNReac = nreac;
  fNEFstat = nfstat;
  fNFstat = fNEbins * fNReac;
  fEmin = TPartIndex::I()->Emin();
  fEmax = TPartIndex::I()->Emax();
  fEGrid = TPartIndex::I()->EGrid();
  fEilDelta = (fNEbins - 1) / log(fEmax / fEmin);
  fFstat = vecfs;

  for (int np = 0; np < nreac; ++np) {
    fRdict[dict[np]] = np;
    fRmap[np] = dict[np];
  }
  // consistency
  for (int i = 0; i < fNReac; ++i)
    if (fRdict[fRmap[i]] != i)
      Fatal("SetPart", "Dictionary mismatch for!");
  return kTRUE;
}

//_________________________________________________________________________
bool TPFstate::SetFinState(int ibin, int reac, const int npart[], const float weight[], const float kerma[],
                           const float en[], const char surv[], const int pid[], const float mom[]) {
  int rnumber = fRdict[reac];
  int ipoint = rnumber * fNEbins + ibin;
  fFstat[ipoint].SetFinState(fNEFstat, npart, weight, kerma, en, surv, pid, mom);
  return kTRUE;
}

//______________________________________________________________________________
bool TPFstate::SampleReac(int preac, float en, int &npart, float &weight, float &kerma, float &enr, const int *&pid,
                          const float *&mom, int &ebinindx) const {
  int rnumber = fRdict[preac];
  if (rnumber == -1) {
    kerma = en;
    npart = 0;
    pid = 0;
    mom = 0;
    ebinindx = -1;
    return kFALSE;
  } else {
    kerma = en;
    double eta = gRandom->Rndm();
    en = en < fEGrid[fNEbins - 1] ? en : fEGrid[fNEbins - 1] * 0.999;
    en = max<double>(en, fEGrid[0]);
    int ibin = log(en / fEGrid[0]) * fEilDelta;
    ibin = ibin < fNEbins - 1 ? ibin : fNEbins - 2;
    double en1 = fEGrid[ibin];
    double en2 = fEGrid[ibin + 1];
    //    if(en1>en || en2<en) {
    //      Error("SetFinState","Wrong bin %d in interpolation: should be %f < %f < %f\n",
    //            ibin, en1, en, en2);
    //      return kFALSE;
    //    }
    double xrat = (en2 - en) / (en2 - en1);
    if (eta > xrat)
      ++ibin;
    ebinindx = ibin;
    int ipoint = rnumber * fNEbins + ibin;
    // in case of any problems with the fstate sampling the primary will be
    // stopped so be prepared for this case and set kerma = en;
    // kerma = en;
    return fFstat[ipoint].SampleReac(npart, weight, kerma, enr, pid, mom);
  }
}

//______________________________________________________________________________
bool TPFstate::SampleReac(int preac, float en, int &npart, float &weight, float &kerma, float &enr, const int *&pid,
                          const float *&mom, int &ebinindx, double randn1, double randn2) const {
  int rnumber = fRdict[preac];
  if (rnumber == -1) {
    kerma = en;
    npart = 0;
    pid = 0;
    mom = 0;
    ebinindx = -1;
    return kFALSE;
  } else {
    // double eta = randn1;
    en = en < fEGrid[fNEbins - 1] ? en : fEGrid[fNEbins - 1] * 0.999;
    en = max<double>(en, fEGrid[0]);
    int ibin = log(en / fEGrid[0]) * fEilDelta;
    ibin = ibin < fNEbins - 1 ? ibin : fNEbins - 2;
    double en1 = fEGrid[ibin];
    double en2 = fEGrid[ibin + 1];
    //    if(en1>en || en2<en) {
    //      Error("SetFinState","Wrong bin %d in interpolation: should be %f < %f < %f\n",
    //            ibin, en1, en, en2);
    //      return kFALSE;
    //    }
    double xrat = (en2 - en) / (en2 - en1);
    if (randn1 > xrat)
      ++ibin;
    int ipoint = rnumber * fNEbins + ibin;
    ebinindx = ibin;
    // in case of any problems with the fstate sampling the primary will be
    // stopped so be prepared for this case and set kerma = en;
    kerma = en;
    return fFstat[ipoint].SampleReac(npart, weight, kerma, enr, pid, mom, randn2);
  }
}

//______________________________________________________________________________
bool TPFstate::SampleRestCaptFstate(int &npart, float &weight, float &kerma, float &enr, const int *&pid,
                                    const float *&mom) const {
  if (fRestCaptFstat) {
    return fRestCaptFstat->SampleReac(npart, weight, kerma, enr, pid, mom);
  } else {
    kerma = 0;
    npart = 0;
    pid = 0;
    mom = 0;
    return kFALSE;
  }
}

//______________________________________________________________________________
bool TPFstate::SampleRestCaptFstate(int &npart, float &weight, float &kerma, float &enr, const int *&pid,
                                    const float *&mom, double randn) const {
  if (fRestCaptFstat) {
    return fRestCaptFstat->SampleReac(npart, weight, kerma, enr, pid, mom, randn);
  } else {
    kerma = 0;
    npart = 0;
    pid = 0;
    mom = 0;
    return kFALSE;
  }
}

//______________________________________________________________________________
bool TPFstate::GetReac(int preac, float en, int ifs, int &npart, float &weight, float &kerma, float &enr,
                       const int *&pid, const float *&mom) const {
  int rnumber = fRdict[preac];
  if (rnumber == -1) {
    kerma = en;
    npart = 0;
    pid = 0;
    mom = 0;
    return kFALSE;
  } else {
    en = en < fEGrid[fNEbins - 1] ? en : fEGrid[fNEbins - 1] * 0.999;
    en = max<double>(en, fEGrid[0]);
    int ibin = log(en / fEGrid[0]) * fEilDelta;
    ibin = ibin < fNEbins - 1 ? ibin : fNEbins - 2;
    double en1 = fEGrid[ibin];
    double en2 = fEGrid[ibin + 1];
    //    if(en1>en || en2<en) {
    //      Error("SetFinState","Wrong bin %d in interpolation: should be %f < %f < %f\n",
    //            ibin, en1, en, en2);
    //      return kFALSE;
    //    }
    if (en - en1 > en2 - en)
      ++ibin;
    int ipoint = rnumber * fNEbins + ibin;
    // in case of any problems with the fstate sampling the primary will be
    // stopped so be prepared for this case and set kerma = en;
    kerma = en;
    return fFstat[ipoint].GetReac(ifs, npart, weight, kerma, enr, pid, mom);
  }
}

#ifdef USE_ROOT
//______________________________________________________________________________
void TPFstate::Streamer(TBuffer &R__b) {
  // Stream an object of class TPFstate.

  if (R__b.IsReading()) {
    R__b.ReadClassBuffer(TPFstate::Class(), this);
    // add the energy grid
    if (!TPartIndex::I()->EGrid()) {
      gFile->Get("PartIndex");
    }
    fEGrid = TPartIndex::I()->EGrid();
  } else {
    R__b.WriteClassBuffer(TPFstate::Class(), this);
  }
}
#endif

//_________________________________________________________________________
void TPFstate::Print(const char *) const {
  printf("Particle=%d Number of x-secs=%d between %g and %g GeV", fPDG, fNReac, fEGrid[0], fEGrid[fNEbins - 1]);
}

//_________________________________________________________________________
bool TPFstate::Resample() {
  if (fVerbose)
    printf("Resampling %s from \nemin = %8.2g emacs = %8.2g, nbins = %d to \n"
           "emin = %8.2g emacs = %8.2g, nbins = %d\n",
           Name(), fEmin, fEmax, fNEbins, TPartIndex::I()->Emin(), TPartIndex::I()->Emax(), TPartIndex::I()->NEbins());
  // Build the original energy grid
  double edelta = exp(1 / fEilDelta);
  double *oGrid = new double[fNEbins];
  double en = fEmin;
  for (int ie = 0; ie < fNEbins; ++ie) {
    oGrid[ie] = en;
    en *= edelta;
  }
  // Build new arrays
  int oNEbins = fNEbins;
  fNEbins = TPartIndex::I()->NEbins();
  // Temporary bins holder
  int *obins = new int[fNEbins];
  fNEFstat = fNEbins * fNReac;
  fEmin = TPartIndex::I()->Emin();
  fEmax = TPartIndex::I()->Emax();
  double oEilDelta = fEilDelta;
  fEilDelta = TPartIndex::I()->EilDelta();
  fEGrid = TPartIndex::I()->EGrid();
  TFinState *oFstat = fFstat;
  fFstat = new TFinState[fNEFstat];

  for (int ien = 0; ien < fNEbins; ++ien) {
    en = fEGrid[ien];
    en = en < oGrid[oNEbins - 1] ? en : oGrid[oNEbins - 1] * 0.999;
    en = max<double>(en, oGrid[0]);
    int ibin = log(fEGrid[ien] / oGrid[0]) * oEilDelta;
    ibin = ibin < oNEbins - 1 ? ibin : oNEbins - 2;
    double en1 = oGrid[ibin];
    double en2 = oGrid[ibin + 1];
    //      if(en1>en || en<en) {
    //	 Error("Interp","Wrong bin %d in interpolation: should be %f < %f < %f\n",
    //	       ibin, en1, en, en2);
    //      }
    double xrat = (en2 - en) / (en2 - en1);
    if (xrat < 0.5)
      obins[ien] = ibin;
    else
      obins[ien] = ibin + 1;
  }

  for (int ir = 0; ir < fNReac; ++ir) {
    int ibase = ir * fNEbins;
    for (int ien = 0; ien < fNEbins; ++ien) {
      fFstat[ibase + ien] = oFstat[obins[ien]];
    }
  }

  delete[] obins;
  delete[] oFstat;
  delete[] oGrid;
  return kTRUE;
}
