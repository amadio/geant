#ifndef TNudySampling_H
#define TNudySampling_H

#define PI acos(-1.0)
class TNudyEndfNuPh;
class TNudyEndfFissionYield;
class TNudyEndfEnergy;
class TNudyEndfEnergyAng;
class TNudyEndfAng;
class TNudyEndfRecoPoint;

class  TNudySampling {

public: 
  TNudySampling ();
  virtual ~TNudySampling ();
private:
  TNudyEndfAng *recoAng;
  TNudyEndfEnergy *recoEnergy;
  TNudyEndfEnergyAng *recoEnergyAng;
  TNudyEndfNuPh *recoNuPh;
  TNudyEndfFissionYield *recoFissY;
  TNudyEndfRecoPoint *recoPoint;

#ifdef USE_ROOT
  ClassDef(TNudySampling, 1) // class for sampling
#endif
};
#endif
