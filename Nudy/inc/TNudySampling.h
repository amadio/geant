#ifndef TNudySampling_H
#define TNudySampling_H

#include "Particle.h"
#define PI acos(-1.0)
#ifdef USE_ROOT
#include "Rtypes.h"
class TRandom;
#endif

class  TNudySampling{

public: 
  TNudySampling ();
  TNudySampling (Particle*, TNudyEndfRecoPoint *recoPoint);
  virtual ~TNudySampling ();
private:
  
 #ifdef USE_ROOT
  TRandom *fRnd;
#endif
 
#ifdef USE_ROOT
  ClassDef(TNudySampling, 1) // class for sampling
#endif

};
#endif
