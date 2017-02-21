#include "TPDecay.h"
#include "TFinState.h"
#include "TPartIndex.h"
#ifndef VECCORE_CUDA
#ifdef USE_ROOT
#include "TBuffer.h"
#endif
#endif

//___________________________________________________________________
TPDecay::TPDecay()
    : fNPart(0), fCTauPerMass(nullptr), fDecay(nullptr), fDecayP(nullptr)
{
}

//___________________________________________________________________
TPDecay::TPDecay(int /*nsample*/, int npart, TFinState *decay) :
   fNPart(npart), fCTauPerMass(nullptr), fDecay(decay), fDecayP(nullptr)
{
   fDecayP = new TFinState*[fNPart];
   for(auto i=0; i<fNPart; ++i)
      fDecayP[i] = &fDecay[i];
}

//___________________________________________________________________
TPDecay::TPDecay(const TPDecay& other) :
   fNPart(other.fNPart), fCTauPerMass(other.fCTauPerMass), fDecay(other.fDecay), fDecayP(other.fDecayP)
{
   fDecayP = new TFinState*[fNPart];
   for(auto i=0; i<fNPart; ++i)
      fDecayP[i] = &fDecay[i];
}

//___________________________________________________________________
TPDecay::~TPDecay() {
   //   delete [] fCTauPerMass;
   //   delete [] fDecay;
   //   if(fDecayP)
   //for(auto i=0; i<fNPart; ++i) delete fDecayP[i];
   delete [] fDecayP;
}

//___________________________________________________________________
bool TPDecay::SampleDecay(int pindex, int &npart, const int *&pid, const float *&mom) const {
  float kerma;
  float weight;
  float en;
  return fDecayP[pindex]->SampleReac(npart, weight, kerma, en, pid, mom);
}

//___________________________________________________________________
bool TPDecay::GetDecay(int pindex, int ifs, int &npart, const int *&pid, const float *&mom) const {
  float kerma;
  float weight;
  float en;
  return fDecayP[pindex]->GetReac(ifs, npart, weight, kerma, en, pid, mom);
}

//___________________________________________________________________
bool TPDecay::HasDecay(int pindex) const {
  if (fDecayP[pindex]->GetNsecs() == 0)
    return false;

  return true;
}

//___________________________________________________________________
void TPDecay::SetCTauPerMass(double *ctaupermass, int np) {
  if (!fCTauPerMass)
    delete fCTauPerMass;
  fCTauPerMass = new double[np];
  for (int ip = 0; ip < np; ++ip)
    fCTauPerMass[ip] = ctaupermass[ip];
}
#ifndef VECCORE_CUDA
#ifdef USE_ROOT
//______________________________________________________________________________
void TPDecay::Streamer(TBuffer &R__b) {
  // Stream an object of class TPFstate.

  if (R__b.IsReading()) {
    R__b.ReadClassBuffer(TPDecay::Class(), this);

    if(fDecayP != nullptr)
      for(auto i=0; i<fNPart; ++i)
        delete fDecayP[i];
    delete [] fDecayP;
    fDecayP =  new TFinState*[fNPart];
    for(auto i=0; i< fNPart; ++i) fDecayP[i]=&fDecay[i];

  } else {
    R__b.WriteClassBuffer(TPDecay::Class(), this);
  }
}
#endif
#endif

//___________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int TPDecay::SizeOf() const {
   size_t size = sizeof(*this);
   size += fNPart * sizeof(double);
   for(auto i=0; i<fNPart; ++i)
      size += fDecayP[i]->SizeOf();
  size -= sizeof(char); // fStore already holds one char
  size = sizeof(double)*((size-1)/sizeof(double)+1);
  return (int) size;
}

//___________________________________________________________________
void TPDecay::Compact() {
   char *start = fStore;
   memcpy(start, fCTauPerMass, fNPart*sizeof(double));
   //   delete [] fCTauPerMass;
   //   fCTauPerMass = (double *) start;
   start += fNPart * sizeof(double);
   for(auto i=0; i<fNPart; ++i) {
      TFinState *px = new(start) TFinState(*fDecayP[i]);
      px->Compact();
      // This line can be reactivated when we remove the back compat array
      //      fFstatP[i]->~TFinState();
      start += px->SizeOf();
      fDecayP[i]=px;
   }
}

//___________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TPDecay::RebuildClass() {
  if(((unsigned long) this) % sizeof(double) != 0) {
    Geant::Fatal("TPDecay::RebuildClass","%s","the class is misaligned\n");
    return;
  }
   char *start = fStore;
   // we consider that the pointer to the final states is stale because it has been read from
   // the file. If this is not the case, this is a leak...
   fCTauPerMass = (double *) start;
   start += fNPart*sizeof(double);
   // we consider that the pointer to the final states is stale because it has been read from
   // the file. If this is not the case, this is a leak...
   fDecayP = new TFinState*[fNPart];
   for(auto i=0; i<fNPart; ++i) {
#ifdef MAGIC_DEBUG
      if(((TFinState*) start)->GetMagic() != -777777) {
        Geant::Fatal("TFinState::RebuildClass","Broken magic 2 %d",((TFinState*) start)->GetMagic());
        return;
      }
#endif
      ((TFinState *) start)->RebuildClass();
      fDecayP[i] = (TFinState *) start;
      if(!fDecayP[i]->CheckAlign()) 
#ifndef VECCORE_CUDA
         exit(1);
#else
         return; 
#endif

      start += ((TFinState*) start)->SizeOf();
   }
   CheckAlign();
}

//___________________________________________________________________
size_t TPDecay::MakeCompactBuffer(char* &b) {
   // First calculate how much we need
   size_t totsize = SizeOf();
   if(b == nullptr) {
      b = (char*) malloc(totsize);
      memset(b,0,totsize);
   }
   TPDecay* dc = new(b) TPDecay(*this);
   dc->Compact();
   return totsize;
}
