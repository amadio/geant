#ifndef TNudyEndfFissionYield_H
#define TNudyEndfFissionYield_H

#include "TNudyEndfRecoPoint.h"
typedef std::vector<double> rowd;
typedef std::vector<rowd > matrixd2;

class  TNudyEndfFissionYield : public TNudyEndfRecoPoint {

public: 
  TNudyEndfFissionYield ();
  TNudyEndfFissionYield (TNudyEndfFile *file);
  virtual ~TNudyEndfFissionYield ();
private:

  double A, AWR, ABN, QX;                                               // standard ENDF parameters
  rowd ein, einc;				        // incident energy
  matrixd2 zafp,fps,zafpc,fpsc,yi, dyi, yc, dyc;// charge, mass, yield (independent and cummulative)
  rowd zafp1,fps1,zafpc1,fpsc1,yi1, dyi1, yc1, dyc1;     // charge, mass, yield (independent and cummulative)
  ClassDef(TNudyEndfFissionYield, 1)                                    // class for an ENDF fission yield reconstruction
};
#endif
