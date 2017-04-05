#include "Geant/Config.h"
#include "Geant/Error.h"
#include "base/Global.h"
#include "base/MessageLogger.h"

#include "TPartIndex.h"
#include "TPFstate.h"
#include "TPartIndex.h"
#include "TEFstate.h"
#include "TPDecay.h"
using vecgeom::kAvogadro;

#ifndef VECCORE_CUDA
#ifdef USE_ROOT
#include "TFile.h"
#endif
#endif


#ifndef VECCORE_CUDA

TEFstate *TEFstate::fElements[NELEM] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

int TEFstate::fNLdElems = 0;
TPDecay* TEFstate::fDecay = nullptr;
#else
TEFstate *fEFElementsDev[NELEM] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

int fEFNLdElemsDev = 0;
TEFstate *fEFElementsHost[NELEM] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

int fEFNLdElemsHost = 0;
TPDecay  *fDecayDev =nullptr;           //! decay table
TPDecay  *fDecayHost=nullptr;           //! decay table
#endif

//___________________________________________________________________
VECCORE_ATT_HOST_DEVICE
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

#ifndef VECCORE_CUDA
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
#endif

//___________________________________________________________________
bool TEFstate::AddPart(int kpart, int pdg, int nfstat, int nreac, const int dict[]) {
  return fPFstateP[kpart]->SetPart(pdg, nfstat, nreac, dict);
}

//___________________________________________________________________
bool TEFstate::AddPart(int kpart, int pdg, int nfstat, int nreac, const int dict[], TFinState vecfs[]) {
  if (!fNEFstat)
    fNEFstat = nfstat;
  else if (fNEFstat != nfstat) {
    Geant::Fatal("TEFstate::AddPart", "Number of final sample states changed during run from %d to %d", fNEFstat, nfstat);
    return 0;
  }
  return fPFstateP[kpart]->SetPart(pdg, nfstat, nreac, dict, vecfs);
}

//___________________________________________________________________
bool TEFstate::AddPartFS(int kpart, int ibin, int reac, const int npart[], const float weight[], const float kerma[],
                         const float en[], const char surv[], const int pid[], const float mom[]) {
  return fPFstateP[kpart]->SetFinState(ibin, reac, npart, weight, kerma, en, surv, pid, mom);
}

//_____________________________________________________________________________
void TEFstate::SetRestCaptFstate(int kpart, const TFinState &fstate) {
   fPFstateP[kpart]->SetRestCaptFstate(fstate);
}

//______________________________________________________________________________
bool TEFstate::HasRestCapture(int partindex) {
  if (partindex < TPartIndex::I()->NPartReac())
    return fPFstateP[partindex]->HasRestCaptFstat();
  return false;
}

//______________________________________________________________________________
bool TEFstate::SampleRestCaptFstate(int kpart, int &npart, float &weight, float &kerma, float &enr, const int *&pid,
                                    const float *&mom) const {
  if (kpart < TPartIndex::I()->NPartReac()) {
    return fPFstateP[kpart]->SampleRestCaptFstate(npart, weight, kerma, enr, pid, mom);
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
    return fPFstateP[kpart]->SampleRestCaptFstate(npart, weight, kerma, enr, pid, mom, randn);
  } else {
    kerma = 0;
    npart = 0;
    pid = 0;
    mom = 0;
    return false;
  }
}

//___________________________________________________________________
VECCORE_ATT_HOST_DEVICE
bool TEFstate::SampleReac(int pindex, int preac, float en, int &npart, float &weight, float &kerma, float &enr,
                          const int *&pid, const float *&mom, int &ebinindx) const {
  return fPFstateP[pindex]->SampleReac(preac, en, npart, weight, kerma, enr, pid, mom, ebinindx);
}

//___________________________________________________________________
VECCORE_ATT_HOST_DEVICE
bool TEFstate::SampleReac(int pindex, int preac, float en, int &npart, float &weight, float &kerma, float &enr,
                          const int *&pid, const float *&mom, int &ebinindx, double randn1, double randn2) const {
  return fPFstateP[pindex]->SampleReac(preac, en, npart, weight, kerma, enr, pid, mom, ebinindx, randn1, randn2);
}

//___________________________________________________________________
  VECCORE_ATT_HOST_DEVICE
bool TEFstate::GetReac(int pindex, int preac, float en, int ifs, int &npart, float &weight, float &kerma, float &enr,
                       const int *&pid, const float *&mom) const {
  return fPFstateP[pindex]->GetReac(preac, en, ifs, npart, weight, kerma, enr, pid, mom);
}

//___________________________________________________________________
#ifdef VECCORE_CUDA 
VECCORE_ATT_HOST_DEVICE
TEFstate *TEFstate::GetElement(int z, int a) {
  int ecode = z * 10000 + a * 10;
#ifdef VECCORE_CUDA_DEVICE_COMPILATION 
  for (int el = 0; el < fEFNLdElemsDev; ++el)
    if (ecode == fEFElementsDev[el]->Ele())
      return fEFElementsDev[el];
  return 0;
#else
  for (int el = 0; el < fEFNLdElemsHost; ++el)
    if (ecode == fEFElementsHost[el]->Ele())
      return fEFElementsHost[el];
  return 0;
#endif  
}
#else   
#ifdef USE_ROOT
TEFstate *TEFstate::GetElement(int z, int a, TFile *f) {
  //   printf("Getting Element %d %d %d\n",z,a,fNLdElems);
  int ecode = z * 10000 + a * 10;
  for (int el = 0; el < fNLdElems; ++el)
    if (ecode == fElements[el]->Ele())
      return fElements[el];

// Element not found in memory, getting it from file
  TFile *ff = gFile;
  if (f)
    ff = f;
  if (!ff)
    Geant::Fatal("TEFstate::GetElement", "No file open!");
  fElements[fNLdElems] = (TEFstate *)ff->Get(TPartIndex::I()->EleSymb(z));
  if (!fElements[fNLdElems]) {
    Geant::Fatal("GetElement", "Element z %d a %d not found", z, a);
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
}
#endif
#endif

//___________________________________________________________________
bool TEFstate::Prune() {
  for (int ip = 0; ip < fNRpart; ++ip)
    fPFstateP[ip]->Prune();
  return true;
}

//___________________________________________________________________
bool TEFstate::Resample() {
  for (int ip = 0; ip < fNRpart; ++ip)
     fPFstateP[ip]->Resample();
  fEmin = TPartIndex::I()->Emin();
  fEmax = TPartIndex::I()->Emax();
  fNEbins = TPartIndex::I()->NEbins();
  fEGrid = TPartIndex::I()->EGrid();
  return true;
}

//___________________________________________________________________
#ifndef VECCORE_CUDA
void TEFstate::Draw(const char * /*option*/) {}
#endif

//___________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int TEFstate::SizeOf() const {
   size_t size = sizeof(*this);
   for(auto i=0; i<fNRpart; ++i)
      size += fPFstateP[i]->SizeOf();
  size -= sizeof(char); // fStore already holds one char
  size = sizeof(double)*((size-1)/sizeof(double)+1);
  return (int) size;
}

//___________________________________________________________________
void TEFstate::Compact() {
   char *start = fStore;
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
VECCORE_ATT_HOST_DEVICE
void TEFstate::RebuildClass() {
   if(((unsigned long) this) % sizeof(double) != 0) {
     Geant::Fatal("TEFstate::RebuildClass","%s","the class is misaligned");
     return;
   }
   char *start = fStore;
   // we consider that the pointer to the final states is stale because it has been read from
   // the file. If this is not the case, this is a leak...
   fPFstateP = new TPFstate*[fNRpart];
   for(auto i=0; i<fNRpart; ++i) {
#ifdef MAGIC_DEBUG
      if(((TPFstate*) start)->GetMagic() != -777777) {
        Geant::Fatal("TEFstate::RebuildClass","Broken magic %d",((TPFstate*) start)->GetMagic());
        return;
      }
#endif
      ((TPFstate *) start)->RebuildClass();
      fPFstateP[i] = (TPFstate *) start;
      if(!fPFstateP[i]->CheckAlign()) 
#ifndef VECCORE_CUDA
        exit(1);
#else
        return;
#endif
      start += ((TPFstate*) start)->SizeOf();
   }
}

//___________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int TEFstate::SizeOfStore() {
   // First calculate how much we need
   int totsize = 0;
#ifndef VECCORE_CUDA
   for(auto i=0; i<fNLdElems; ++i) totsize += fElements[i]->SizeOf();
#else
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
   for(auto i=0; i<fEFNLdElemsDev; ++i) totsize += fEFElementsDev[i]->SizeOf();
#else
   for(auto i=0; i<fEFNLdElemsHost; ++i) totsize += fEFElementsHost[i]->SizeOf();
#endif
#endif
   return totsize + 2 * sizeof(int);
}

//___________________________________________________________________
int TEFstate::MakeCompactBuffer(char* &b) {
   // First calculate how much we need
   size_t totsize = SizeOfStore();
   if(b == nullptr) {
      // Now allocate buffer
      b = (char*) malloc(totsize);
      memset(b,0,totsize);
   }
   char* start = b;
   memcpy(start,&totsize,sizeof(int));
   start += sizeof(int);
#ifndef  VECCORE_CUDA
   memcpy(start,&fNLdElems,sizeof(int));
   start += sizeof(int);
   // now copy and compact
   for(auto i=0; i<fNLdElems; ++i) {
      TEFstate *el = new(start) TEFstate(*fElements[i]);
      el->Compact();
      start += el->SizeOf();
   }
 #endif
   return totsize;
}

#ifndef VECCORE_CUDA
#define EXIT_OR_RETURN(val) exit(val);
#else
#define EXIT_OR_RETURN(val) return;
#endif

//___________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TEFstate::RebuildStore(char *b) {
  typedef TEFstate *StateArray_t[NELEM];
#ifndef VECCORE_CUDA
  int &nElems(fNLdElems);
  StateArray_t &elements(fElements);
#else
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
  int &nElems(fEFNLdElemsDev);
  StateArray_t &elements(fEFElementsDev);
#else
  int &nElems(fEFNLdElemsHost);
  StateArray_t &elements(fEFElementsHost);
#endif
#endif
  char *start = b;
  long size = 0;
  memcpy(&size,start,sizeof(int));
  start += sizeof(int);
  memcpy(&nElems,start,sizeof(int));
  start += sizeof(int);
  for(auto i=0; i<nElems; ++i) {
    TEFstate *current = (TEFstate *) start;
#ifdef MAGIC_DEBUG
    if(current->GetMagic() != -777777) {
      Geant::Fatal("TEFstate::RebuildStore","Broken Magic %d\n",current->GetMagic());
      EXIT_OR_RETURN(1);
    }
#endif
    current->RebuildClass();
    elements[i] = current;
    if(!elements[i]->CheckAlign())
      EXIT_OR_RETURN(1);
    start += current->SizeOf();
  }
  if((start - b) != size) {
    Geant::Fatal("TEFstate::RebuildStore","expected size %ld  found size %ld\n",size,start - b);
    EXIT_OR_RETURN(1);
  }
}
