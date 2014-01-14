#include "GeantTrackStat.h"
#include "GeantTrack.h"

ClassImp(GeantTrackStat)

//______________________________________________________________________________
GeantTrackStat::GeantTrackStat(Int_t nslots)
      :TObject(), fNslots(nslots), fNtracks(0), fNsteps(0), fMutex()
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
void GeantTrackStat::AddTrack(const GeantTrack &track)
{
// Add a track
   fNtracks[track.fEvslot]++;
   fNsteps[track.fEvslot] += track.fNsteps;
}   

//______________________________________________________________________________
void GeantTrackStat::AddTrack(const GeantTrack_v &trackv, Int_t itr)
{
// Add a track from an array
   fNtracks[trackv.fEvslotV[itr]]++;
   fNsteps[trackv.fEvslotV[itr]] += trackv.fNstepsV[itr];
}   

//______________________________________________________________________________
void GeantTrackStat::AddTracks(const GeantTrack_v &trackv)
{
// Remove statistics for tracks
   fMutex.Lock();
   Int_t ntracks = trackv.GetNtracks();
   for (Int_t i=0; i<ntracks; i++) {
      fNtracks[trackv.fEvslotV[i]]++;
      fNsteps[trackv.fEvslotV[i]] += trackv.fNstepsV[i];
   }   
   fMutex.UnLock();
}   

//______________________________________________________________________________
void GeantTrackStat::RemoveTracks(const GeantTrack_v &trackv)
{
// Remove statistics for tracks
   fMutex.Lock();
   Int_t ntracks = trackv.GetNtracks();
   for (Int_t i=0; i<ntracks; i++) {
      // do *NOT* remove new tracks since they were not added yet anywhere
      if (trackv.fStatusV[i]==kNew) continue;
      fNtracks[trackv.fEvslotV[i]]--;
      fNsteps[trackv.fEvslotV[i]] -= trackv.fNstepsV[i];
   }   
   fMutex.UnLock();
}   

//______________________________________________________________________________
GeantTrackStat & GeantTrackStat::operator+=(const GeantTrackStat &other)
{
// Compound addition.
   if (fNslots != other.fNslots) {
      Error("operator+=","Different number of slots");
      return *this;
   }
   for (Int_t i=0; i<fNslots; i++) {
      fNtracks[i] += other.fNtracks[i];
      fNsteps[i] += other.fNsteps[i];
   }
   return *this;
}      
   
//______________________________________________________________________________
GeantTrackStat & GeantTrackStat::operator-=(const GeantTrackStat &other)
{
// Compound addition.
   if (fNslots != other.fNslots) {
      Error("operator+=","Different number of slots");
      return *this;
   }
   for (Int_t i=0; i<fNslots; i++) {
      fNtracks[i] -= other.fNtracks[i];
      fNsteps[i] -= other.fNsteps[i];
   }
   return *this;
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
void GeantTrackStat::Print(Option_t *) const
{
// Print statistics
   printf("slot: ");
   for (Int_t i=0; i<fNslots; i++) {
      printf("%5d ", i);
   }
   printf("\n");
   printf("ntr:  ");
   for (Int_t i=0; i<fNslots; i++) {
      printf("%5d ", fNtracks[i]);
   }
   printf("\n");
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
