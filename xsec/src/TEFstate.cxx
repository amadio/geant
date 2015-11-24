#ifdef USE_ROOT
#include "TFile.h"
#include "TError.h"
#endif

#include <TPFstate.h>
#include <TPartIndex.h>
#include "TEFstate.h"
#include "base/Global.h"
#include "base/MessageLogger.h"
using vecgeom::kAvogadro;

TEFstate *TEFstate::fElements[NELEM] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

int TEFstate::fNLdElems = 0;

//___________________________________________________________________
TEFstate::TEFstate() :
   fEGrid(TPartIndex::I()->EGrid()),
   fAtcm3(0), 
   fEmin(0), 
   fEmax(0), 
   fEilDelta(0), 
   fDens(0), 
   fEle(0), 
   fNEbins(0), 
   fNEFstat(0), 
   fNRpart(0), 
   fPFstate(nullptr),
   fPFstateP(nullptr)
{}

//___________________________________________________________________
TEFstate::TEFstate(int z, int a, float dens) :
   fEGrid(TPartIndex::I()->EGrid()), 
   fAtcm3(dens * kAvogadro * 1e-24 / TPartIndex::I()->WEle(z)),
   fEmin(TPartIndex::I()->Emin()), 
   fEmax(TPartIndex::I()->Emax()), 
   fEilDelta(TPartIndex::I()->EilDelta()), 
   fDens(dens), 
   fEle(z * 10000 + a * 10), 
   fNEbins(TPartIndex::I()->NEbins()),
   fNEFstat(0),
   fNRpart(TPartIndex::I()->NPartReac()), 
   fPFstate(new TPFstate[fNRpart]),
   fPFstateP(new TPFstate*[fNRpart])
{
   for(auto i=0; i< fNRpart; ++i) fPFstateP[i] = &fPFstate[i];
}

//___________________________________________________________________
TEFstate::TEFstate(const TEFstate &other) :
   fEGrid(TPartIndex::I()->EGrid()), 
   fAtcm3(other.fAtcm3),
   fEmin(other.fEmin),
   fEmax(other.fEmax),
   fEilDelta(other.fEilDelta),
   fDens(other.fDens),
   fEle(other.fEle),
   fNEbins(other.fNEbins),
   fNEFstat(other.fNEFstat),
   fNRpart(other.fNRpart),
   fPFstate(other.fPFstate),
   fPFstateP(other.fPFstateP)
{
}

//___________________________________________________________________
TEFstate::~TEFstate() 
{ 
   if(fPFstateP != nullptr) 
      for(auto i=0; i<fNRpart; ++i) delete fPFstateP[i];
   delete [] fPFstateP;
   delete [] fPFstate; 
}

#ifdef USE_ROOT
//______________________________________________________________________________
void TEFstate::Streamer(TBuffer &R__b)
{
   // Stream an object of class TEXsec.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TEFstate::Class(),this);
      if(fPFstateP != nullptr) 
	 for(auto ipart=0; ipart<fNRpart; ++ipart) delete fPFstateP[ipart];
      delete [] fPFstateP; 
      fPFstateP = new TPFstate*[fNRpart];
      for(auto i=0; i<fNRpart; ++i) fPFstateP[i] = &fPFstate[i];
   } else {
      R__b.WriteClassBuffer(TEFstate::Class(),this);
   }
}
#endif

//___________________________________________________________________
bool TEFstate::AddPart(int kpart, int pdg, int nfstat, int nreac, const int dict[]) {
  return fPFstate[kpart].SetPart(pdg, nfstat, nreac, dict);
}

//___________________________________________________________________
bool TEFstate::AddPart(int kpart, int pdg, int nfstat, int nreac, const int dict[], TFinState vecfs[]) {
  if (!fNEFstat)
    fNEFstat = nfstat;
  else if (fNEFstat != nfstat) {
    log_fatal(std::cout, "AddPart", "Number of final sample states changed during run from %d to %d", fNEFstat, nfstat);
    exit(1);
  }
  return fPFstate[kpart].SetPart(pdg, nfstat, nreac, dict, vecfs);
}

//___________________________________________________________________
bool TEFstate::AddPartFS(int kpart, int ibin, int reac, const int npart[], const float weight[], const float kerma[],
                         const float en[], const char surv[], const int pid[], const float mom[]) {
  return fPFstate[kpart].SetFinState(ibin, reac, npart, weight, kerma, en, surv, pid, mom);
}

//_____________________________________________________________________________
void TEFstate::SetRestCaptFstate(int kpart, const TFinState &fstate) { fPFstate[kpart].SetRestCaptFstate(fstate); }

//______________________________________________________________________________
bool TEFstate::HasRestCapture(int partindex) {
  if (partindex < TPartIndex::I()->NPartReac())
    return fPFstate[partindex].HasRestCaptFstat();
  return false;
}

//______________________________________________________________________________
bool TEFstate::SampleRestCaptFstate(int kpart, int &npart, float &weight, float &kerma, float &enr, const int *&pid,
                                    const float *&mom) const {
  if (kpart < TPartIndex::I()->NPartReac()) {
    return fPFstate[kpart].SampleRestCaptFstate(npart, weight, kerma, enr, pid, mom);
  } else {
    kerma = 0;
    npart = 0;
    pid = 0;
    mom = 0;
    return false;
  }
}

//______________________________________________________________________________
bool TEFstate::SampleRestCaptFstate(int kpart, int &npart, float &weight, float &kerma, float &enr, const int *&pid,
                                    const float *&mom, double randn) const {
  if (kpart < TPartIndex::I()->NPartReac()) {
    return fPFstate[kpart].SampleRestCaptFstate(npart, weight, kerma, enr, pid, mom, randn);
  } else {
    kerma = 0;
    npart = 0;
    pid = 0;
    mom = 0;
    return false;
  }
}

//___________________________________________________________________
bool TEFstate::SampleReac(int pindex, int preac, float en, int &npart, float &weight, float &kerma, float &enr,
                          const int *&pid, const float *&mom, int &ebinindx) const {
  return fPFstate[pindex].SampleReac(preac, en, npart, weight, kerma, enr, pid, mom, ebinindx);
}

//___________________________________________________________________
bool TEFstate::SampleReac(int pindex, int preac, float en, int &npart, float &weight, float &kerma, float &enr,
                          const int *&pid, const float *&mom, int &ebinindx, double randn1, double randn2) const {
  return fPFstate[pindex].SampleReac(preac, en, npart, weight, kerma, enr, pid, mom, ebinindx, randn1, randn2);
}

//___________________________________________________________________
bool TEFstate::GetReac(int pindex, int preac, float en, int ifs, int &npart, float &weight, float &kerma, float &enr,
                       const int *&pid, const float *&mom) const {
  return fPFstate[pindex].GetReac(preac, en, ifs, npart, weight, kerma, enr, pid, mom);
}

//___________________________________________________________________
TEFstate *TEFstate::GetElement(int z, int a, TFile *f) {
  //   printf("Getting Element %d %d %d\n",z,a,fNLdElems);
  int ecode = z * 10000 + a * 10;
  for (int el = 0; el < fNLdElems; ++el)
    if (ecode == fElements[el]->Ele())
      return fElements[el];

// Element not found in memory, getting it from file
#ifdef USE_ROOT
  TFile *ff = gFile;
  if (f)
    ff = f;
  if (!ff)
    ::Fatal("TEFstate::GetElement", "No file open!");
  fElements[fNLdElems] = (TEFstate *)ff->Get(TPartIndex::I()->EleSymb(z));
  if (!fElements[fNLdElems]) {
    ::Fatal("GetElement", "Element z %d a %d not found", z, a);
    return 0; // just to make the compiler happy
  } else {
    // We loaded the element, however we have to see whether
    // the energy grid is the right one
    // NO, don't need to check. It will be loaded from xsec.root
    //      if(FloatDiff(TPartIndex::I()->Emin(),fElements[fNLdElems]->Emin(),1e-7) ||
    //	 FloatDiff(TPartIndex::I()->Emax(),fElements[fNLdElems]->Emax(),1e-7) ||
    //	 TPartIndex::I()->NEbins() != fElements[fNLdElems]->NEbins())
    // we have to resize the energy grid of the element
    //	 fElements[fNLdElems]->Resample();
    return fElements[fNLdElems++];
  }
#else
  log_fatal(std::cout, "Element z %d a %d not found\n", z, a);
  exit(1);
  return 0;
#endif
}

//___________________________________________________________________
bool TEFstate::Prune() {
  for (int ip = 0; ip < fNRpart; ++ip)
    fPFstate[ip].Prune();
  return true;
}

//___________________________________________________________________
bool TEFstate::Resample() {
  for (int ip = 0; ip < fNRpart; ++ip)
     fPFstate[ip].Resample();
  fEmin = TPartIndex::I()->Emin();
  fEmax = TPartIndex::I()->Emax();
  fNEbins = TPartIndex::I()->NEbins();
  fEGrid = TPartIndex::I()->EGrid();
  return true;
}

//___________________________________________________________________
void TEFstate::Draw(const char * /*option*/) {}

//___________________________________________________________________
int TEFstate::SizeOf() const {
   size_t size = sizeof(*this);
   for(auto i=0; i<fNRpart; ++i)
      size += fPFstateP[i]->SizeOf();
   return (int) size - sizeof(TPFstate); // fStore already holds one TPXsec
}

//___________________________________________________________________
void TEFstate::Compact() {
   char *start = (char*) fStore;
   for(auto i=0; i<fNRpart; ++i) {
      TPFstate *px = new(start) TPFstate(*fPFstateP[i]);
      px->Compact();
      // This line can be reactivated when we remove the back compat array
      //      fPFstateP[i]->~TPXsec();
      start += px->SizeOf();
      fPFstateP[i]=px;
   }
}

//___________________________________________________________________
void TEFstate::RebuildClass() {
   char *start = (char*) fStore;
   for(auto i=0; i<fNRpart; ++i) {
#ifdef MAGIC_DEBUG
      if(((TPFstate*) start)->GetMagic() != -777777) {
	 cout << "TPFstate::Broken magic " << ((TPFstate*) start)->GetMagic() << endl;
	 exit(1);
      }
#endif
      ((TPFstate *) start)->RebuildClass();
      fPFstateP[i] = (TPFstate *) start;
      start += ((TPFstate*) start)->SizeOf();
   }
}

//___________________________________________________________________
size_t TEFstate::MakeCompactBuffer(char* &b) {
   // First calculate how much we need
   size_t totsize = 0;
   for(auto i=0; i<fNLdElems; ++i) totsize += fElements[i]->SizeOf();
   // Now allocate buffer
   b = (char*) malloc(totsize);
   char* start = b;
   // now copy and compact
   for(auto i=0; i<fNLdElems; ++i) {
      TEFstate *el = new(start) TEFstate(*fElements[i]);
      el->Compact();
      start += el->SizeOf();
   }
   return totsize;
}

//___________________________________________________________________
void TEFstate::RebuildStore(size_t size, int nelem, char *b) {
   fNLdElems = 0;
   char *start = b;
   for(auto i=0; i<nelem; ++i) {
      TEFstate *current = (TEFstate *) start;
#ifdef MAGIC_DEBUG
      if(current->GetMagic() != -777777) {
	 cout << "TEFstate::Broken Magic " << current->GetMagic() << endl;
	 exit(1);
      }
#endif
      current->RebuildClass();
      fElements[fNLdElems++] = current;
      start += current->SizeOf();     
   }
   if((size_t)(start - b) != size) {
      cout << "TEFstate::RebuildStore: expected size " << size 
	   << " found size " << start - b << endl;
      exit(1);
   }
}

