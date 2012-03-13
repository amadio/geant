#include "GeantParticleBuffer.h"

#include "GeantVolumeBasket.h"

ClassImp(GeantParticleBuffer)

//______________________________________________________________________________
GeantParticleBuffer::GeantParticleBuffer(Int_t nthreads, Int_t max_pending)
                    :TObject(),
                     fNthreads(nthreads),
                     fMaxPending(max_pending),
                     fNpending(0),
                     fTrackInd(0),
                     fBaskets(0)
{
// Constructor
   fNpending = new Int_t[nthreads];
   fTrackInd = new Int_t*[nthreads];
   fBaskets = new GeantVolumeBasket**[nthreads];
   for (Int_t tid=0; tid<nthreads; tid++) {
      fNpending[tid] = 0;
      fTrackInd[tid] = new Int_t[max_pending];
      fBaskets[tid] = new GeantVolumeBasket*[max_pending];
      for (Int_t j=0; j<max_pending; j++) {
         fTrackInd[tid][j] = 0;
         fBaskets[tid][j]  = 0;
      }   
   }
}

//______________________________________________________________________________
GeantParticleBuffer::~GeantParticleBuffer()
{
// Destructor.
   for (Int_t i=0; i<fNthreads; i++) {
      delete [] fTrackInd[i];
      delete [] fBaskets[i];
   }
   delete [] fNpending;
   delete [] fTrackInd;
   delete [] fBaskets;
}

//______________________________________________________________________________
void GeantParticleBuffer::AddPendingTrack(Int_t itrack, GeantVolumeBasket *basket, Int_t tid)
{
// Add a pending track by a given thread.
   fTrackInd[tid][fNpending[tid]] = itrack;
   fBaskets[tid][fNpending[tid]] = basket;
   fNpending[tid]++;
}

//______________________________________________________________________________
void GeantParticleBuffer::FlushBaskets()
{
// Flush all baskets. To be called by a single thread.
   for (Int_t tid=0; tid<fNthreads; tid++) {
      for (Int_t itrack=0; itrack<fNpending[tid]; itrack++) {
         fBaskets[tid][itrack]->AddTrack(fTrackInd[tid][itrack]);
      }
      fNpending[tid] = 0;
   }
}   

//______________________________________________________________________________
Bool_t GeantParticleBuffer::IsEmpty() const
{
// Check if buffer is empty.
   Int_t nparticles = 0;
   for (Int_t tid=0; tid<fNthreads; tid++) {
      nparticles += fNpending[tid];
   }
   return (nparticles==0);
}
