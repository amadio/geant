#include "TPXsec.h"
#ifdef USE_ROOT
#include "TFile.h"
#include "TRandom.h"
#endif
#ifdef USE_VECGEOM_NAVIGATOR
#include "base/RNG.h"
using vecgeom::RNG;
#endif

#include "Geant/Error.h"

using std::max;

int TPXsec::fVerbose = 0;

//_________________________________________________________________________
TPXsec::TPXsec()
   : fPDG(0), fNEbins(0), fNCbins(0), fNXsec(0), fNTotXs(0), fNXSecs(0),
     fEGrid(TPartIndex::I()->EGrid()), fMSangle(nullptr), fMSansig(nullptr), fMSlength(nullptr),
     fMSlensig(nullptr), fdEdx(nullptr), fTotXs(nullptr), fXSecs(nullptr),  fEmin(0), fEmax(0), fEilDelta(0) {
   int np = TPartIndex::I()->NProc();
   while (np--)
      fRdict[np] = fRmap[np] = -1;
}

//_________________________________________________________________________
TPXsec::TPXsec(int pdg, int nxsec)
   :  fPDG(pdg), fNEbins(TPartIndex::I()->NEbins()),
      fNCbins(0), fNXsec(nxsec), fNTotXs(fNEbins), fNXSecs(fNEbins * fNXsec),
      fEGrid(TPartIndex::I()->EGrid()), fMSangle(nullptr),
      fMSansig(nullptr), fMSlength(nullptr), fMSlensig(nullptr), fdEdx(nullptr), fTotXs(nullptr),
      fXSecs(nullptr), fEmin(TPartIndex::I()->Emin()), fEmax(TPartIndex::I()->Emax()),
      fEilDelta((TPartIndex::I()->NEbins() - 1) / log(fEmax / fEmin)) {
   int np = TPartIndex::I()->NProc();
   while (np--)
      fRdict[np] = fRmap[np] = -1;
}

//_________________________________________________________________________
TPXsec::TPXsec(const TPXsec &other): fPDG(other.fPDG), fNEbins(other.fNEbins),
				     fNCbins(other.fNCbins), fNXsec(other.fNXsec),
				     fNTotXs(other.fNTotXs), fNXSecs(other.fNXSecs),
				     fEGrid(other.fEGrid),
				     fMSangle(other.fMSangle), fMSansig(other.fMSansig),
				     fMSlength(other.fMSlength), fMSlensig(other.fMSlensig),
				     fdEdx(other.fdEdx), fTotXs(other.fTotXs), fXSecs(other.fXSecs),
				     fEmin(other.fEmin), fEmax(other.fEmax),
				     fEilDelta(other.fEilDelta)
{
   memcpy(fRdict, other.fRdict, FNPROC*sizeof(int));
   memcpy(fRmap, other.fRmap, FNPROC*sizeof(int));
}

//_________________________________________________________________________
TPXsec::~TPXsec() {
  delete[] fMSangle;
  delete[] fMSansig;
  delete[] fMSlength;
  delete[] fMSlensig;
  delete[] fdEdx;
  delete[] fTotXs;
  delete[] fXSecs;
}

//_________________________________________________________________________
int TPXsec::SizeOf() const {
   size_t size = sizeof(*this);
   size += 5 * fNCbins * sizeof(float);
   size += fNTotXs * sizeof(float);
   size += fNXSecs * sizeof(float);
   size -= sizeof(char); // fStore already takes one char
   size = sizeof(double)*((size-1)/sizeof(double)+1);
   return (int) size;
}

//_________________________________________________________________________
void TPXsec::Compact() {
   int size = 0;
   char *start = fStore;
   if(fMSangle) {
      size = fNCbins * sizeof(float);
      memcpy(start, fMSangle, size);
      fMSangle = (float*) start;
      start += size;
   }
   if(fMSansig) {
      size = fNCbins * sizeof(float);
      memcpy(start, fMSansig, size);
      fMSansig = (float*) start;
      start += size;
   }
   if(fMSlength) {
      size = fNCbins * sizeof(float);
      memcpy(start, fMSlength, size);
      fMSlength = (float*) start;
      start += size;
   }
   if(fMSlensig) {
      size = fNCbins * sizeof(float);
      memcpy(start, fMSlensig, size);
      fMSlensig = (float*) start;
      start += size;
   }
   if(fdEdx) {
      size = fNCbins * sizeof(float);
      memcpy(start, fdEdx, size);
      fdEdx = (float*) start;
      start += size;
   }
   if(fTotXs) {
      size = fNTotXs * sizeof(float);
      memcpy(start, fTotXs, size);
      fTotXs = (float*) start;
      start += size;
   }
   if(fXSecs) {
      size = fNXSecs * sizeof(float);
      memcpy(start, fXSecs, size);
      fXSecs = (float*) start;
      start += size;
   }
}

//______________________________________________________________________________
void TPXsec::RebuildClass() {
  if(((unsigned long) this) % sizeof(double) != 0) {
      cout << "TPXsec::RebuildClass: the class is misaligned" << endl;
      exit(1);
  }
   // Reset fEgrid, may be in a different place
   fEGrid = TPartIndex::I()->EGrid();
   int size = 0;
   char *start = fStore;
   if(fMSangle) {
      //      cout << "Original fMSangle " << fMSangle << " new pointer " << start << endl;
      fMSangle = (float*) start;
      size = fNCbins * sizeof(float);
      start += size;
   }
   if(fMSansig) {
      //      cout << "Original fMSansig " << fMSansig << " new pointer " << start << endl;
      fMSansig = (float*) start;
      size = fNCbins * sizeof(float);
      start += size;
   }
   if(fMSlength) {
      //      cout << "Original fMSlength " << fMSlength << " new pointer " << start << endl;
      fMSlength = (float*) start;
      size = fNCbins * sizeof(float);
      start += size;
   }
   if(fMSlensig) {
      //      cout << "Original fMSlensig " << fMSlensig << " new pointer " << start << endl;
      fMSlensig = (float*) start;
      size = fNCbins * sizeof(float);
      start += size;
   }
   if(fdEdx) {
      //      cout << "Original fdEdx " << fdEdx << " new pointer " << start << endl;
      fdEdx = (float*) start;
      size = fNCbins * sizeof(float);
      start += size;
   }
   if(fTotXs) {
      //      cout << "Original fTotXs " << fTotXs << " new pointer " << start << endl;
      fTotXs = (float*) start;
      size = fNTotXs * sizeof(float);
      start += size;
   }
   if(fXSecs) {
      //      cout << "Original fXSecs " << fXSecs << " new pointer " << start << endl;
      fXSecs = (float*) start;
      size = fNXSecs * sizeof(float);
      start += size;
   }
}

#ifdef USE_ROOT
//______________________________________________________________________________
void TPXsec::Streamer(TBuffer &R__b) {
  // Stream an object of class TPXsec.

  if (R__b.IsReading()) {
    R__b.ReadClassBuffer(TPXsec::Class(), this);
    // add the energy grid
    if (!TPartIndex::I()->EGrid()) {
      gFile->Get("PartIndex");
    }
    fEGrid = TPartIndex::I()->EGrid();
  } else {
    R__b.WriteClassBuffer(TPXsec::Class(), this);
  }
}
#endif

//_________________________________________________________________________
void TPXsec::Interp(double egrid[], float value[], int nbins, double eildelta, int stride, double en, float result[]) {
  en = en < egrid[nbins - 1] ? en : egrid[nbins - 1] * 0.999;
  en = max<double>(en, egrid[0]);
  int ibin = log(en / egrid[0]) * eildelta;
  ibin = ibin < nbins - 1 ? ibin : nbins - 2;
  double en1 = egrid[ibin];
  double en2 = egrid[ibin + 1];
  //   if(en1>en || en2<en) {
  //      Error("Interp","Wrong bin %d in interpolation: should be %f < %f < %f\n",
  //	    ibin, en1, en, en2);
  //   }
  double xrat = (en2 - en) / (en2 - en1);
  for (int ival = 0; ival < stride; ++ival) {
    result[ival] = xrat * value[ibin * stride + ival] + (1 - xrat) * value[(ibin + 1) * stride + ival];
  }
}

//___________________________________________________________________
bool TPXsec::Prune() {
  //
  // Dangerous business... delete unwanted cross-sections
  //
  delete[] fTotXs;
  fTotXs = 0;
  fNTotXs = 0;
  delete[] fMSangle;
  fMSangle = 0;
  delete[] fMSansig;
  fMSansig = 0;
  delete[] fMSlength;
  fMSlength = 0;
  delete[] fMSlensig;
  fMSlensig = 0;
  delete[] fdEdx;
  fdEdx = 0;
  fNCbins = 0;
  return true;
}

//___________________________________________________________________
bool TPXsec::Resample() {
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
  fNTotXs = fNEbins = TPartIndex::I()->NEbins();
  if (fNCbins)
    fNCbins = fNEbins;
  fNXSecs = fNEbins * fNXsec;
  fEmin = TPartIndex::I()->Emin();
  fEmax = TPartIndex::I()->Emax();
  double oEilDelta = fEilDelta;
  fEilDelta = TPartIndex::I()->EilDelta();
  fEGrid = TPartIndex::I()->EGrid();
  float *lXSecs = new float[fNXSecs];
  float *lTotXs = new float[fNEbins];
  // Total x-secs and partial channels
  for (int ien = 0; ien < fNEbins; ++ien) {
    Interp(oGrid, fXSecs, oNEbins, oEilDelta, fNXsec, fEGrid[ien], &lXSecs[ien * fNXsec]);
    double xnorm = 0;
    // recheck normalisation
    for (int ixs = 0; ixs < fNXsec; ++ixs)
      xnorm += lXSecs[ien * fNXsec + ixs];
    xnorm = 1 / xnorm;
    for (int ixs = 0; ixs < fNXsec; ++ixs)
      lXSecs[ien * fNXsec + ixs] *= xnorm;
    Interp(oGrid, fTotXs, oNEbins, oEilDelta, 1, fEGrid[ien], &lTotXs[ien]);
  }
  delete[] fXSecs;
  fXSecs = lXSecs;
  delete[] fTotXs;
  fTotXs = lTotXs;
  // Only for charged particles
  if (fNCbins) {
    float *lMSangle = new float[fNCbins];
    float *lMSansig = new float[fNCbins];
    float *lMSlength = new float[fNCbins];
    float *lMSlensig = new float[fNCbins];
    float *ldEdx = new float[fNCbins];
    for (int ien = 0; ien < fNEbins; ++ien) {
      Interp(oGrid, fMSangle, oNEbins, oEilDelta, 1, fEGrid[ien], &lMSangle[ien]);
      Interp(oGrid, fMSansig, oNEbins, oEilDelta, 1, fEGrid[ien], &lMSansig[ien]);
      Interp(oGrid, fMSlength, oNEbins, oEilDelta, 1, fEGrid[ien], &lMSlength[ien]);
      Interp(oGrid, fMSlensig, oNEbins, oEilDelta, 1, fEGrid[ien], &lMSlensig[ien]);
      Interp(oGrid, fdEdx, oNEbins, oEilDelta, 1, fEGrid[ien], &ldEdx[ien]);
    }
    delete[] fMSangle;
    fMSangle = lMSangle;
    delete[] fMSansig;
    fMSansig = lMSansig;
    delete[] fMSlength;
    fMSlength = lMSlength;
    delete[] fMSlensig;
    fMSlensig = lMSlensig;
    delete[] fdEdx;
    fdEdx = ldEdx;
  }
  delete[] oGrid;
  return true;
}

//___________________________________________________________________
bool TPXsec::SetPart(int pdg, int nxsec) {
  fPDG = pdg;
  fNEbins = TPartIndex::I()->NEbins();
  fNXsec = nxsec;
  fNTotXs = fNEbins;
  fNXSecs = fNEbins * fNXsec;
  fEmin = TPartIndex::I()->Emin();
  fEmax = TPartIndex::I()->Emax();
  fEGrid = TPartIndex::I()->EGrid();
  fEilDelta = TPartIndex::I()->EilDelta();
  return true;
}

//___________________________________________________________________
bool TPXsec::SetPartXS(const float xsec[], const int dict[]) {
  delete[] fXSecs;
  fXSecs = new float[fNXSecs];
  for (int jxsec = 0; jxsec < fNXsec; ++jxsec)
    for (int jbin = 0; jbin < fNEbins; ++jbin)
      fXSecs[jbin * fNXsec + jxsec] = xsec[jxsec * fNEbins + jbin];
  for (int i = 0; i < TPartIndex::I()->NProc(); ++i)
    fRdict[i] = fRmap[i] = -1;
  for (int i = 0; i < fNXsec; ++i) {
    fRdict[dict[i]] = i;
    fRmap[i] = dict[i];
  }
  // consistency
  for (int i = 0; i < fNXsec; ++i)
    if (fRdict[fRmap[i]] != i)
      Geant::Fatal("TPXsec::SetPartXS", "Dictionary mismatch for!");

  delete[] fTotXs;
  fTotXs = new float[fNTotXs];
  for (int i = 0; i < fNEbins; ++i) {
    fTotXs[i] = 0;
    for (int j = 0; j < fNXsec; ++j)
      fTotXs[i] += fXSecs[i * fNXsec + j];
    if (fTotXs[i])
      for (int j = 0; j < fNXsec; ++j)
        fXSecs[i * fNXsec + j] /= fTotXs[i];
  }
  return true;
}

//___________________________________________________________________
bool TPXsec::SetPartIon(const float dedx[]) {
  delete[] fdEdx;
  fNCbins = fNEbins;
  fdEdx = new float[fNCbins];
  memcpy(fdEdx, dedx, fNCbins * sizeof(float));
  return true;
}

//___________________________________________________________________
bool TPXsec::SetPartMS(const float angle[], const float ansig[], const float length[], const float lensig[]) {
  fNCbins = fNEbins;

  delete[] fMSangle;
  fMSangle = new float[fNCbins];
  memcpy(fMSangle, angle, fNCbins * sizeof(float));

  delete[] fMSansig;
  fMSansig = new float[fNCbins];
  memcpy(fMSansig, ansig, fNCbins * sizeof(float));

  delete[] fMSlength;
  fMSlength = new float[fNCbins];
  memcpy(fMSlength, length, fNCbins * sizeof(float));

  delete[] fMSlensig;
  fMSlensig = new float[fNCbins];
  memcpy(fMSlensig, lensig, fNCbins * sizeof(float));

  return true;
}

//_________________________________________________________________________
void TPXsec::Print(const char *) const {
  printf("Particle=%d Number of x-secs=%d between %g and %g GeV", fPDG, fNXsec, fEGrid[0], fEGrid[fNEbins - 1]);
}

//_________________________________________________________________________
float TPXsec::DEdx(double en) const {
  if (!fdEdx)
    return 0;
  en = en < fEGrid[fNEbins - 1] ? en : fEGrid[fNEbins - 1] * 0.999;
  en = max<double>(en, fEGrid[0]);
  int ibin = log(en / fEGrid[0]) * fEilDelta;
  ibin = ibin < fNEbins - 1 ? ibin : fNEbins - 2;
  //   double en1 = fEmin*exp(ibin/fEilDelta);
  //   double en2 = fEmin*exp((ibin+1)/fEilDelta);
  double en1 = fEGrid[ibin];
  double en2 = fEGrid[ibin + 1];
  //   if(en1>en || en2<en) {
  //      Error("XS","Wrong bin %d in interpolation: should be %f < %f < %f\n",
  //	    ibin, en1, en, en2);
  //      return 0;
  //   }
  double xrat = (en2 - en) / (en2 - en1);
  double dedx = xrat * fdEdx[ibin] + (1 - xrat) * fdEdx[ibin + 1];
  /*   printf("ibin %d en1 %f en %f en2 %f xs1 %f xs2 %f xrat %f xsec %f\n",
       ibin,en1,en,en2,xs1,xs2,xrat,xsec); */
  return dedx;
}

//_________________________________________________________________________
bool TPXsec::MS(double en, float &ang, float &asig, float &len, float &lsig) const {
  if (!fNCbins) {
    ang = asig = len = lsig = 0;
    return false;
  }
  en = en < fEGrid[fNEbins - 1] ? en : fEGrid[fNEbins - 1] * 0.999;
  en = max<double>(en, fEGrid[0]);
  int ibin = log(en / fEGrid[0]) * fEilDelta;
  ibin = ibin < fNEbins - 1 ? ibin : fNEbins - 2;
  //   double en1 = fEmin*exp(ibin/fEilDelta);
  //  double en2 = fEmin*exp((ibin+1)/fEilDelta);
  double en1 = fEGrid[ibin];
  double en2 = fEGrid[ibin + 1];
  //   if(en1>en || en2<en) {
  //      Error("XS","Wrong bin %d in interpolation: should be %f < %f < %f\n",
  //	    ibin, en1, en, en2);
  //      return 0;
  //   }
  double xrat = (en2 - en) / (en2 - en1);
  ang = xrat * fMSangle[ibin] + (1 - xrat) * fMSangle[ibin + 1];
  asig = xrat * fMSansig[ibin] + (1 - xrat) * fMSansig[ibin + 1];
  len = xrat * fMSlength[ibin] + (1 - xrat) * fMSlength[ibin + 1];
  lsig = xrat * fMSlensig[ibin] + (1 - xrat) * fMSlensig[ibin + 1];
  return true;
}

//_________________________________________________________________________
int TPXsec::SampleReac(double en) const {
  en = en < fEGrid[fNEbins - 1] ? en : fEGrid[fNEbins - 1] * 0.999;
  en = max<double>(en, fEGrid[0]);
  int ibin = log(en / fEGrid[0]) * fEilDelta;
  ibin = ibin < fNEbins - 1 ? ibin : fNEbins - 2;
  //   double en1 = fEmin*exp(ibin/fEilDelta);
  // double en2 = fEmin*exp((ibin+1)/fEilDelta);
  double en1 = fEGrid[ibin];
  double en2 = fEGrid[ibin + 1];
  //   if(en1>en || en2<en) {
  //      Error("SampleReac","Wrong bin %d in interpolation: should be %f < %f < %f\n",
  //	    ibin, en1, en, en2);
  //      return -1;
  //   }

  double xrat = (en2 - en) / (en2 - en1);
  double xnorm = 1.;
  while (1) {
#ifdef USE_ROOT
    double ran = xnorm * gRandom->Rndm();
#else
    double ran = RNG::Instance().uniform();
#endif
    double xsum = 0;
    for (int i = 0; i < fNXsec; ++i) {
      xsum += xrat * fXSecs[ibin * fNXsec + i] + (1 - xrat) * fXSecs[(ibin + 1) * fNXsec + i];
      if (ran <= xsum)
        return fRmap[i];
    }
    xnorm = xsum;
  }
}

//_________________________________________________________________________
int TPXsec::SampleReac(double en, double randn) const {
  en = en < fEGrid[fNEbins - 1] ? en : fEGrid[fNEbins - 1] * 0.999;
  en = max<double>(en, fEGrid[0]);
  int ibin = log(en / fEGrid[0]) * fEilDelta;
  ibin = ibin < fNEbins - 1 ? ibin : fNEbins - 2;
  //   double en1 = fEmin*exp(ibin/fEilDelta);
  // double en2 = fEmin*exp((ibin+1)/fEilDelta);
  double en1 = fEGrid[ibin];
  double en2 = fEGrid[ibin + 1];
  //   if(en1>en || en2<en) {
  //      Error("SampleReac","Wrong bin %d in interpolation: should be %f < %f < %f\n",
  //	    ibin, en1, en, en2);
  //      return -1;
  //   }

  double xrat = (en2 - en) / (en2 - en1);
  double xnorm = 1.;
  while (1) {
    double ran = xnorm * randn;
    double xsum = 0;
    for (int i = 0; i < fNXsec; ++i) {
      xsum += xrat * fXSecs[ibin * fNXsec + i] + (1 - xrat) * fXSecs[(ibin + 1) * fNXsec + i];
      if (ran <= xsum)
        return fRmap[i];
    }
    xnorm = xsum;
  }
}

//_________________________________________________________________________
bool TPXsec::XS_v(int npart, int rindex, const double en[], double lam[]) const {
  //   printf("fEGrid %p\n",fEGrid);
  double ene;
  for (int ip = 0; ip < npart; ++ip) {
    ene = en[ip] < fEGrid[fNEbins - 1] ? en[ip] : fEGrid[fNEbins - 1] * 0.999;
    ene = max<double>(en[ip], fEGrid[0]);
    int ibin = log(ene / fEGrid[0]) * fEilDelta;
    ibin = ibin < fNEbins - 1 ? ibin : fNEbins - 2;
    //   double en1 = fEmin*exp(ibin/fEilDelta);
    //   double en2 = fEmin*exp((ibin+1)/fEilDelta);
    double en1 = fEGrid[ibin];
    double en2 = fEGrid[ibin + 1];
    //    if(en1>ene || en2<ene) {
    //      Error("XS","Wrong bin %d in interpolation: should be %f < %f < %f\n",
    //            ibin, en1, ene, en2);
    //      return 0;
    //    }
    double xrat = (en2 - ene) / (en2 - en1);

    double xtot = (xrat * fTotXs[ibin] + (1 - xrat) * fTotXs[ibin + 1]);
    double xsec = 1;
    if (rindex < TPartIndex::I()->NProc() - 1) {
      int rnumber = fRdict[rindex];
      if (rnumber < 0) {
        Geant::Error("TPXsec::XS", "No %s for %s\n", TPartIndex::I()->ProcName(rindex),
                     TPartIndex::I()->PartName(TPartIndex::I()->PartIndex(fPDG)));
        return -1;
      }
      xsec = xrat * fXSecs[ibin * fNXsec + rnumber] + (1 - xrat) * fXSecs[(ibin + 1) * fNXsec + rnumber];
    }
    lam[ip] = xsec * xtot;
  }
  return true;
}

//_________________________________________________________________________
float TPXsec::XS(int rindex, double en, bool verbose) const {
  //   printf("fEGrid %p\n",fEGrid);
  en = en < fEGrid[fNEbins - 1] ? en : fEGrid[fNEbins - 1] * 0.999;
  en = max<double>(en, fEGrid[0]);
  int ibin = log(en / fEGrid[0]) * fEilDelta;
  ibin = ibin < fNEbins - 1 ? ibin : fNEbins - 2;
  //   double en1 = fEmin*exp(ibin/fEilDelta);
  //   double en2 = fEmin*exp((ibin+1)/fEilDelta);
  double en1 = fEGrid[ibin];
  double en2 = fEGrid[ibin + 1];
  //   if(en1>en || en2<en) {
  //      Error("XS","Wrong bin %d in interpolation: should be %f < %f < %f\n",
  //	    ibin, en1, en, en2);
  //      return 0;
  //   }
  double xrat = (en2 - en) / (en2 - en1);

  double xtot = (xrat * fTotXs[ibin] + (1 - xrat) * fTotXs[ibin + 1]);
  double xsec = 1;
  if (rindex < TPartIndex::I()->NProc() - 1) {
    int rnumber = fRdict[rindex];
    if (rnumber < 0 && verbose) {
      Geant::Error("TPXsec::XS", "No %s for %s\n", TPartIndex::I()->ProcName(rindex),
                   TPartIndex::I()->PartName(TPartIndex::I()->PartIndex(fPDG)));
      return -1;
    }
    xsec = xrat * fXSecs[ibin * fNXsec + rnumber] + (1 - xrat) * fXSecs[(ibin + 1) * fNXsec + rnumber];
  }
  return xsec * xtot;
}

//_________________________________________________________________________
void TPXsec::Dump() const {
  printf("Particle %d NXsec %d emin %f emax %f NEbins %d ElDelta %f\n", fPDG, fNXsec, fEGrid[0], fEGrid[fNEbins - 1],
         fNEbins, fEilDelta);
  printf("MSangle %p, MSlength %p, dEdx %p, TotXs %p, XSecs %p, Rdict %p\n", (void *)fMSangle, (void *)fMSlength,
         (void *)fdEdx, (void *)fTotXs, (void *)fXSecs, (const void *)fRdict);
}
