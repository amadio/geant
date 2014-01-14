#ifndef GEANT_TRACK_STAT
#define GEANT_TRACK_STAT

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifndef ROOT_TMutex
#include "TMutex.h"
#endif

//______________________________________________________________________________
// Statistics object for a track array
//______________________________________________________________________________
class GeantTrack;
class GeantTrack_v;

//______________________________________________________________________________
class GeantTrackStat : public TObject
{
public:
   Int_t    fNslots;    // Number of event slots
   Int_t   *fNtracks;   //[fNslots] Number of tracks from an event
   Int_t   *fNsteps;    //[fNslots] Cumulated number of steps per event
   TMutex   fMutex;     //! Mutex
   
private:
   GeantTrackStat(const GeantTrackStat &other);
   GeantTrackStat &operator=(const GeantTrackStat &other);

public:
   GeantTrackStat() : TObject(), fNslots(0), fNtracks(0), fNsteps(0), fMutex() {}
   GeantTrackStat(Int_t nslots);
   virtual ~GeantTrackStat();
   
   // operators +=, -=
   GeantTrackStat &operator+=(const GeantTrackStat &other);
   GeantTrackStat &operator-=(const GeantTrackStat &other);
   
   void       AddTrack(const GeantTrack &track);
   void       AddTrack(const GeantTrack_v &trackv, Int_t itr);
   void       AddTracks(const GeantTrack_v &trackv);
   void       RemoveTracks(const GeantTrack_v &trackv);
   void       InitArrays(Int_t nslots);
   void       Print(Option_t *option="") const;
   void       Reset();
   
   ClassDef(GeantTrackStat, 1)      // Track statistics
};
#endif
