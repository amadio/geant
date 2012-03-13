#ifndef GEANT_VOLUMEBASKET
#define GEANT_VOLUMEBASKET

#ifndef ROOT_TObject
#include "TObject.h"
#endif
#include "GeantTrack.h"
 
#ifndef ROOT_TGeoVolume
#include "TGeoVolume.h"
#endif

class TGeoBranchArray;

//______________________________________________________________________________
class GeantVolumeBasket : public TObject {
protected:
   TGeoVolume       *fVolume;                // Volume for which applies
   Int_t             fNtracks;               // Number of tracks
   Int_t             fFirstFree;             // First un-processed track
   Int_t             fMaxTracks;             // Max number of tracks
   Int_t            *fIndex;                 //[fNtracks] Track indices in the global stack

public:
   GeantVolumeBasket(TGeoVolume *vol);
   virtual ~GeantVolumeBasket();
   
   void              AddTrack(Int_t itrack);
   virtual void      Clear(Option_t *option="");
   void              ComputeTransportLengthSingle(Int_t *trackin);
   void              ComputeTransportLength(Int_t ntracks, Int_t *trackin);
   TGeoBranchArray  *GetBranchArray(GeantTrack *track) const {return track->path;}
   TGeoBranchArray  *GetBranchArray(Int_t itrack) const;
   Int_t             GetIndex(Int_t itrack) const {return fIndex[itrack];}
   const Int_t      *GetIndArray() const          {return fIndex;}
   Int_t             GetNchunks(Int_t nthreads) const;
   Int_t             GetNtotal() const            {return fNtracks;}
   const char       *GetName() const              {return (fVolume)?fVolume->GetName():ClassName();}
   Int_t             GetNtracks() const           {return fNtracks;}
   GeantTrack       *GetTrack(Int_t itrack) const;
   TGeoVolume       *GetVolume() const            {return fVolume;}
   void              GetWorkload(Int_t &indmin, Int_t &indmax);
   virtual void      Print(Option_t *option="") const;
   Bool_t            PropagateTrack(Int_t *trackin);
   void              PropagateTracks(Int_t ntracks, Int_t *trackin, Int_t &nout, Int_t *trackout, Int_t &ntodo, Int_t *tracktodo, Int_t &ncross, Int_t *trackcross);
   static void       ResetStep(Int_t ntracks, Int_t *array);
   void              TransportSingle();
//   void              TransportTracks();
   
   ClassDef(GeantVolumeBasket,1)  // A path in geometry represented by the array of indices
};   
#endif
