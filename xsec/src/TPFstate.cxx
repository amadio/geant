#include "TPFstate.h"
#include "TFinState.h"
#include "Geant/Error.h"
#ifdef USE_ROOT
#include <TMath.h>
#include <TFile.h>
#include <TRandom.h>
#else
#include "base/RNG.h"
using vecgeom::RNG;
#endif
using std::max;

int TPFstate::fVerbose = 0;

#include "Geant/Error.h"

//_________________________________________________________________________
TPFstate::TPFstate()
  : fNEbins(0),
    fNEFstat(0),
    fNFstat(0),
    fNReac(0),
    fFstat(nullptr),
    fFstatP(nullptr),
    fRestCaptFstat(nullptr),
    fEGrid(TPartIndex::I()->EGrid()),
    fEmin(0),
    fEmax(0),
    fEilDelta(0),
    fPDG(0)
{
  int np = TPartIndex::I()->NProc();
  while (np--)
    fRdict[np] = fRmap[np] = -1;
}

//_________________________________________________________________________
TPFstate::TPFstate(int pdg, int nfstat, int nreac, const int dict[]) :
   fNEbins(TPartIndex::I()->NEbins()),
   fNEFstat(nfstat),
   fNFstat(fNEbins * fNEFstat),
   fNReac(nreac),
   fFstat(new TFinState[fNFstat]),
   fFstatP(new TFinState*[fNFstat]),
   fRestCaptFstat(nullptr),
   fEGrid(TPartIndex::I()->EGrid()),
   fEmin(TPartIndex::I()->Emin()),
   fEmax(TPartIndex::I()->Emax()),
   fEilDelta((fNEbins - 1) / log(fEmax / fEmin)),
   fPDG(pdg)
{
   for(auto i=0; i< fNFstat; ++i) fFstatP[i]=&fFstat[i];
   int np = TPartIndex::I()->NProc();
   while (np--) {
      fRdict[dict[np]] = np;
      fRmap[np] = dict[np];
   }
   // consistency
   for (int i = 0; i < fNReac; ++i)
      if (fRdict[fRmap[i]] != i)
	 Geant::Fatal("TPFstate::TPFstate", "Dictionary mismatch for!");
}

//_________________________________________________________________________
TPFstate::TPFstate(const TPFstate& other) :
   fNEbins(other.fNEbins),
   fNEFstat(other.fNEFstat),
   fNFstat(other.fNFstat),
   fNReac(other.fNReac),
   fFstat(other.fFstat),
   fFstatP(other.fFstatP),
   fRestCaptFstat(other.fRestCaptFstat),
   fEGrid(TPartIndex::I()->EGrid()),
   fEmin(other.fEmin),
   fEmax(other.fEmax),
   fEilDelta(other.fEilDelta),
   fPDG(other.fPDG)
{
  int np = TPartIndex::I()->NProc();
  memcpy(fRdict,other.fRdict,np*sizeof(int));
  memcpy(fRmap,other.fRmap,np*sizeof(int));
}

//_________________________________________________________________________
TPFstate::~TPFstate() {
  delete [] fFstatP;
  delete[] fFstat;
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
      Geant::Fatal("TPFstate::SetPart", "Dictionary mismatch for!");
  return true;
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
      Geant::Fatal("TPFstateSetPart", "Dictionary mismatch for!");
  return true;
}

//_________________________________________________________________________
bool TPFstate::SetFinState(int ibin, int reac, const int npart[], const float weight[], const float kerma[],
                           const float en[], const char surv[], const int pid[], const float mom[]) {
  int rnumber = fRdict[reac];
  int ipoint = rnumber * fNEbins + ibin;
  fFstatP[ipoint]->SetFinState(fNEFstat, npart, weight, kerma, en, surv, pid, mom);
  return true;
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
    return false;
  } else {
    kerma = en;
#ifdef USE_ROOT
    double eta = gRandom->Rndm();
#else
    double eta = RNG::Instance().uniform();
#endif
    en = en < fEGrid[fNEbins - 1] ? en : fEGrid[fNEbins - 1] * 0.999;
    en = max<double>(en, fEGrid[0]);
    int ibin = log(en / fEGrid[0]) * fEilDelta;
    ibin = ibin < fNEbins - 1 ? ibin : fNEbins - 2;
    double en1 = fEGrid[ibin];
    double en2 = fEGrid[ibin + 1];
    //    if(en1>en || en2<en) {
    //      Error("SetFinState","Wrong bin %d in interpolation: should be %f < %f < %f\n",
    //            ibin, en1, en, en2);
    //      return false;
    //    }
    double xrat = (en2 - en) / (en2 - en1);
    if (eta > xrat)
      ++ibin;
    ebinindx = ibin;
    int ipoint = rnumber * fNEbins + ibin;
    // in case of any problems with the fstate sampling the primary will be
    // stopped so be prepared for this case and set kerma = en;
    // kerma = en;
    return fFstatP[ipoint]->SampleReac(npart, weight, kerma, enr, pid, mom);
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
    return false;
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
    //      return false;
    //    }
    double xrat = (en2 - en) / (en2 - en1);
    if (randn1 > xrat)
      ++ibin;
    int ipoint = rnumber * fNEbins + ibin;
    ebinindx = ibin;
    // in case of any problems with the fstate sampling the primary will be
    // stopped so be prepared for this case and set kerma = en;
    kerma = en;
    return fFstatP[ipoint]->SampleReac(npart, weight, kerma, enr, pid, mom, randn2);
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
    return false;
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
    return false;
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
    return false;
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
    //      return false;
    //    }
    if (en - en1 > en2 - en)
      ++ibin;
    int ipoint = rnumber * fNEbins + ibin;
    // in case of any problems with the fstate sampling the primary will be
    // stopped so be prepared for this case and set kerma = en;
    kerma = en;
    return fFstatP[ipoint]->GetReac(ifs, npart, weight, kerma, enr, pid, mom);
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

    if(fFstatP != nullptr)
      for(auto i=0; i<fNFstat; ++i)
        delete fFstatP[i];
    delete [] fFstatP;
    fFstatP =  new TFinState*[fNFstat];
    for(auto i=0; i< fNFstat; ++i) fFstatP[i]=&fFstat[i];

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
  return true;
}

//___________________________________________________________________
int TPFstate::SizeOf() const {
   size_t size = sizeof(*this);
   if(fRestCaptFstat != nullptr) size += fRestCaptFstat->SizeOf();
   for(auto i=0; i<fNFstat; ++i)
      size += fFstatP[i]->SizeOf();
    size -= sizeof(char); // fStore already holds one TPXsec
    size = sizeof(double)*((size-1)/sizeof(double)+1);
    return (int) size;
}

//___________________________________________________________________
void TPFstate::Compact() {
   char *start = fStore;
   if(fRestCaptFstat != nullptr) {
      TFinState *px = new(start) TFinState(*fRestCaptFstat);
      px->Compact();
      fRestCaptFstat->~TFinState();
      start += px->SizeOf();
      fRestCaptFstat=px;
   }

   for(auto i=0; i<fNFstat; ++i) {
      TFinState *px = new(start) TFinState(*fFstatP[i]);
      px->Compact();
      // This line can be reactivated when we remove the back compat array
      //      fFstatP[i]->~TFinState();
      start += px->SizeOf();
      fFstatP[i]=px;
   }
}

//___________________________________________________________________
void TPFstate::RebuildClass() {
  if(((unsigned long) this) % sizeof(double) != 0) {
      cout << "TPFstate::RebuildClass: the class is misaligned" << endl;
      exit(1);
  }
   char *start = fStore;
   // we consider that the pointer energy grid is stale
   fEGrid = TPartIndex::I()->EGrid();
   // we consider that the pointer to the final states is stale because it has been read from
   // the file. If this is not the case, this is a leak...
   fFstatP = new TFinState*[fNFstat];
   if(fRestCaptFstat != nullptr) {
#ifdef MAGIC_DEBUG
      if(((TFinState*) start)->GetMagic() != -777777) {
	 cout << "TFinState::Broken magic 1 " << ((TFinState*) start)->GetMagic() << endl;
	 exit(1);
      }
#endif
      ((TFinState *) start)->RebuildClass();
      fRestCaptFstat = (TFinState *) start;
      if(!fRestCaptFstat->CheckAlign()) exit(1);
      start += ((TFinState*) start)->SizeOf();
   }
   for(auto i=0; i<fNFstat; ++i) {
#ifdef MAGIC_DEBUG
      if(((TFinState*) start)->GetMagic() != -777777) {
	 cout << "TFinState::Broken magic 2 " << ((TFinState*) start)->GetMagic() << endl;
	 exit(1);
      }
#endif
      ((TFinState *) start)->RebuildClass();
      fFstatP[i] = (TFinState *) start;
      fFstatP[i]->CheckAlign();
      start += ((TFinState*) start)->SizeOf();
   }
}
