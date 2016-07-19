#ifndef TNudyEndfEnergyAng_H
#define TNudyEndfEnergyAng_H

#include "TNudyEndfRecoPoint.h"

#define PI acos(-1.0)
typedef std::vector<double> rowd;
typedef std::vector<int> rowint;
typedef std::vector<rowd > matrixd2;
typedef std::vector<std::vector<rowd > > matrixd3;
typedef std::vector<std::vector<std::vector<rowd > > > matrixd4;
#ifdef USE_ROOT
#include "Rtypes.h"
class TRandom3;
#endif

class  TNudyEndfEnergyAng : public TNudyEndfRecoPoint {

public: 
  TNudyEndfEnergyAng ();
  TNudyEndfEnergyAng (TNudyEndfFile *file, double []);
  virtual double GetCos64(int elemid, int mt, double energyK);
  virtual double GetCos6(int elemid, int mt, double energyK);
  virtual double GetEnergy6(int elemid, int mt, double energyK);
  virtual int GetLaw6(int ielemId, int mt);
  virtual int GetZd6(int ielemId, int mt);
  virtual int GetAd6(int ielemId, int mt);
  virtual ~TNudyEndfEnergyAng ();
private:
  double recursionLinearFile5Prob(double x1, double x2, double pdf1, double pdf2);
  double recursionLinearFile4(int i, double x1, double x2, double pdf1, double pdf2);
  double recursionLinearProb(double x1, double x2, double pdf1, double pdf2);
  void Law2OnlyAngle();

  double A, AWR, ABN, QX;                         // standard ENDF parameters
  int testv;
  int NR, NP;                         // standard ENDF parameters for range and interpolation
  int NR2, NE;                            // standard ENDF parameters for range and interpolation
  int NR3, NE2;				// standard ENDF parameters for range and interpolation
//  double *INorm, QValue[999];				// ENDF parameter and Q values from file 3
  rowd cosFile4, energyFile5;
  rowd cosPdfFile4, energyPdfFile5;
  rowd cosCdfFile4, energyCdfFile5;
  rowd edes6,f06,r6,a6;
  rowd fE1, fP1, fE2, fP2, fE3, fP3, INorm;
  rowint zd, ad;
  rowint law;				// law6 numbers
  rowint MtNumbers,MtNumbers6, MtNumbers4;				// MT numbers
  rowint MtLct;				// LCT numbers
  rowint nbt1,int1;
  int nr1, np1;                         // standard ENDF parameters
  rowint nbt2,int2;
  int nr2, np2;                         // standard ENDF parameters
  rowint nbt3,int3;
  int nr3, np3;                         // standard ENDF parameters
  rowd ein,cosc,cdfc,pdfc,lCoef1,cosin,cosinpdf,cosincdf;
  rowd eoute,cdfe,pdfe;
  matrixint Mt6Values;             // MT values 
  matrixd2 cos2d,cosinpdf2d,cosincdf2d,cos2dc,pdf2dc,cdf2dc,lCoef,ein2d,ein2dc;
  matrixd3 cos3d,cosinpdf3d,cosincdf3d,cos3dc,pdf3dc,cdf3dc;
  matrixd2 eout2de,pdf2de,cdf2de;
  matrixd3 eout3de,pdf3de,cdf3de;
  matrixd4 eout4de,pdf4de,cdf4de;
#ifdef USE_ROOT
  TRandom3 *fRnd;
#endif
#ifdef USE_ROOT
  ClassDef(TNudyEndfEnergyAng, 1) // class for an ENDF reconstruction
#endif
};
#endif
