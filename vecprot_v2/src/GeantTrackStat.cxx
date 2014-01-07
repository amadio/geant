#include "GeantTrackStat.h"

ClassImp(GeantTrackStat)

//______________________________________________________________________________
GeantTrackStat::GeantTrackStat(Int_t nslots)
      :TObject(), fNslots(nslots), fNtracks(0), fNsteps(0)
{
// Ctor
   InitArrays(nslots);
}

//______________________________________________________________________________
GeantTrackStat::~GeantTrackStat()
{
// Dtor.
   delete [] fNtracks;
   delete [] fNsteps;
}
   
//______________________________________________________________________________
void GeantTrackStat::InitArrays(Int_t nslots)
{
// Initialize arrays
   fNslots = nslots;
   delete [] fNtracks;
   delete [] fNsteps;
   fNtracks = new Int_t[nslots];
   memset(fNtracks, 0, nslots*sizeof(Int_t));
   fNsteps  = new Int_t[nslots];
   memset(fNsteps, 0, nslots*sizeof(Int_t));
}

//______________________________________________________________________________
void GeantTrackStat::Reset()
{
// Reset statistics
   if (fNslots) {
      memset(fNtracks, 0, fNslots*sizeof(Int_t));
      memset(fNsteps, 0, fNslots*sizeof(Int_t));
   }
}  
