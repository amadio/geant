#ifndef GEANT_MYHIT
#define GEANT_MYHIT

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

//______________________________________________________________________________
class MyHit {
public:
   double fX;      // position
   double fY;
   double fZ;
   double fDe;     // energy loss
   int    fVolId;  // volume Id
   int    fDetId;  // replica (segmentation)
   
   MyHit() : fX(0), fY(0), fZ(0), fDe(0), fVolId(0), fDetId(0) {}
   MyHit(double x, double y, double z, double de, int volid, int detid);
   ~MyHit() {}
   
   void               Reset() {fX=0.; fY=0.; fZ=0.; fDe=0.; fVolId=0; fDetId=0;}
   void               AddHit();
   
   ClassDefNV(MyHit, 1)      // User hit
};
#endif
