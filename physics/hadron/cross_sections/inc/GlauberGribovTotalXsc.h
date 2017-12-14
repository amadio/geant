#ifndef GlauberGribovTotalXsc_H
#define GlauberGribovTotalXsc_H

#include "HadronicCrossSection.h"
#include "Parameterizations.h"

namespace geantphysics {

  class GlauberGribovTotalXsc : public HadronicCrossSection
{
public:

  GlauberGribovTotalXsc ();
  virtual ~GlauberGribovTotalXsc ();
   
  double GetIsotopeCrossSection(int particlePDG, double mass, double energyKin, int Z, int N);
  
 private:
 
  inline double GetParticleBarCorTot(int particlePDG, int Z);
  
  static const double fNeutronBarCorrectionTot[93];
  
  static const double fProtonBarCorrectionTot[93];

  static const double fPionPlusBarCorrectionTot[93];

  static const double fPionMinusBarCorrectionTot[93];

};

/////////////////////////////////////////////////////////////////////
//
// return correction at Tkin = 90*GeV GG -> Barashenkov tot xsc, when it 
// is available, else return 1.0


inline double GlauberGribovTotalXsc::GetParticleBarCorTot(int particlePDG, int Z)
{
  if(Z >= 2 && Z <= 92)
  {
    if(      particlePDG == 2212 ) return fProtonBarCorrectionTot[Z]; 
    else if( particlePDG == 2112 ) return fNeutronBarCorrectionTot[Z]; 
    else if( particlePDG == 211  ) return fPionPlusBarCorrectionTot[Z];
    else if( particlePDG == -211 ) return fPionMinusBarCorrectionTot[Z];
    else return 1.0;
  }
  else return 1.0;
}

}       // namespace geantphysics

#endif  // GlauberGribovTotalXsc_H


// proton = 2212, neutron = 2112
// anti-proton = -2212
// pi+ = 211
// pi- = -211
// K+ = 321, K0L = 130
// K- = -321, K0S = 310
// sigma- = 3112, 
