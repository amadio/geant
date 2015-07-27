#ifndef GEANT_VOLUMEBASKET
#define GEANT_VOLUMEBASKET

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifndef ROOT_TGeoVolume
#include "TGeoVolume.h"
#endif

#ifndef GEANT_BASKET
#include "GeantBasket.h"
#endif

#include "GeantTrack.h"

class TGeoBranchArray;

//______________________________________________________________________________
class GeantVolumeBasket : public TObject {
protected:
  TGeoVolume *fVolume; // Volume for which applies
  int fNumber;         // Number assigned

public:
  GeantVolumeBasket(TGeoVolume *vol, int number);
  virtual ~GeantVolumeBasket();

  const char *GetName() const { return (fVolume) ? fVolume->GetName() : ClassName(); }
  int GetNumber() const { return fNumber; }
  TGeoVolume *GetVolume() const { return fVolume; }
  virtual void Print(const char *option = "") const;

  void ComputeTransportLength(int ntracks, int *trackin);
  void PropagateTracks(int ntracks, int *trackin, int &nout, int *trackout, int &ntodo, int *tracktodo, int &ncross,
                       int *trackcross);
  static void ResetStep(int ntracks, int *array);

  ClassDef(GeantVolumeBasket, 1) // A path in geometry represented by the array of indices
};

#endif // GEANT_VOLUMEBASKET
