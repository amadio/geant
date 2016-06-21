#ifndef TNudyEndfFissionYield_H
#define TNudyEndfFissionYield_H

#include "TNudyEndfRecoPoint.h"

class  TNudyEndfFissionYield : public TNudyEndfRecoPoint {

public: 
  TNudyEndfFissionYield ();
  TNudyEndfFissionYield (TNudyEndfFile *file);
  virtual ~TNudyEndfFissionYield ();
private:

  double A, AWR, ABN, QX;                                               // standard ENDF parameters
  std::vector<double> ein, einc;				        // incident energy
  std::vector<std::vector<double> >zafp,fps,zafpc,fpsc,yi, dyi, yc, dyc;// charge, mass, yield (independent and cummulative)
  std::vector<double> zafp1,fps1,zafpc1,fpsc1,yi1, dyi1, yc1, dyc1;     // charge, mass, yield (independent and cummulative)
  ClassDef(TNudyEndfFissionYield, 1)                                    // class for an ENDF fission yield reconstruction
};
#endif
