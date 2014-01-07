#ifndef GEANT_TRACK_STAT
#define GEANT_TRACK_STAT

#ifndef ROOT_TObject
#include "TObject.h"
#endif

//______________________________________________________________________________
// Statistics object for a track array
//______________________________________________________________________________

//______________________________________________________________________________
class GeantTrackStat : public TObject
{
public:
   Int_t    fNslots;    // Number of event slots
   Int_t   *fNtracks;   //[fNslots] Number of tracks from an event
   Int_t   *fNsteps;    //[fNslots] Cumulated number of steps per event
   
private:
   GeantTrackStat(const GeantTrackStat &other);
   GeantTrackStat &operator=(const GeantTrackStat &other);

public:
   GeantTrackStat() : TObject(), fNslots(0), fNtracks(0), fNsteps(0) {}
   GeantTrackStat(Int_t nslots);
   virtual ~GeantTrackStat();
   
   void       InitArrays(Int_t nslots);
   void       Reset();
   
   ClassDef(GeantTrackStat, 1)      // Track statistics
};
#endif
