#include "MyHit.h"
#include "GeantFactoryStore.h"
#include "TStopwatch.h"

ClassImp(MyHit)

    //______________________________________________________________________________
    MyHit::MyHit(double x, double y, double z, double edep, double time, int event, int track, int volid, int detid)
    : fX(x), fY(y), fZ(z), fEdep(edep), 
      fTime(time), fEvent(event), fTrack(track), fVolId(volid), fDetId(detid) {
  // Ctor..
}

//______________________________________________________________________________
void MyHit::AddHit() {}
