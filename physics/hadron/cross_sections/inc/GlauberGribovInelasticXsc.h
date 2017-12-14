#ifndef GlauberGribovInelasticXsc_H
#define GlauberGribovInelasticXsc_H

#include "HadronicCrossSection.h"
#include "Parameterizations.h"

namespace geantphysics {

class GlauberGribovInelasticXsc : public HadronicCrossSection
{
public:

  GlauberGribovInelasticXsc ();
  virtual ~GlauberGribovInelasticXsc ();
   
  double GetIsotopeCrossSection(int particlePDG, double mass, double energyKin, int Z, int A);
  
 private:
  
  inline double GetParticleBarCorIn(int particlePDG, int Z);
  
  static const double fNeutronBarCorrectionIn[93];
  
  static const double fProtonBarCorrectionIn[93];

  static const double fPionPlusBarCorrectionIn[93];

  static const double fPionMinusBarCorrectionIn[93];

};

/////////////////////////////////////////////////////////////////////
//
// return correction at Tkin = 90*GeV GG -> Barashenkov in xsc, when it 
// is available, else return 1.0


inline double GlauberGribovInelasticXsc::GetParticleBarCorIn(int particlePDG, int Z)
{
  if(Z >= 2 && Z <= 92)
  {
    if(      particlePDG == 2212 ) return fProtonBarCorrectionIn[Z]; 
    else if( particlePDG == 2112 ) return fNeutronBarCorrectionIn[Z]; 
    else if( particlePDG == 211  ) return fPionPlusBarCorrectionIn[Z];
    else if( particlePDG == -211 ) return fPionMinusBarCorrectionIn[Z];
    else return 1.0;
  }
  else return 1.0;
}

}       // namespace geantphysics

#endif  // GlauberGribovInelasticXsc_H


// proton = 2212, neutron = 2112
// anti-proton = -2212
// pi+ = 211
// pi- = -211
// K+ = 321, K0L = 130
// K- = -321, K0S = 310
// sigma- = 3112, 
