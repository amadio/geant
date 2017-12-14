#ifndef HadronNucleusXsc_H
#define HadronNucleusXsc_H

namespace geantphysics {

class HadronNucleusXsc 
{
public:

  HadronNucleusXsc ();
  virtual ~HadronNucleusXsc ();
   
  double CalculateCrossSection(int particlePDG, double mass, double energyKin, int Z, int A);
  
  inline double GetTotalXsc() {return fNucleusTotalXsc;};
  inline double GetElasticXsc() {return fNucleusElasticXsc;};
  inline double GetInelasticXsc() {return fNucleusInelasticXsc;};

 private:
  
  double GetHadronNucleonXscNS(int particlePDG, double mass, double energyKin, int targetPDG);
  double GetHadronNucleonXscPDG(int particlePDG, double mass, double energyKin, int targetPDG);
  double GetKaonNucleonXscGG(int particlePDG, double mass, double energyKin, int targetPDG);

  double CalcMandelstamS(const double mp, const double mt, const double Plab);
  double GetCoulombBarrier(int particlePDG, double proj_mass, double energyKin, int targetPDG, double target_mass);
  double GetNucleusRadius(int At);

  inline double GetParticleBarCorTot(int particlePDG, int Z);
  inline double GetParticleBarCorIn(int particlePDG, int Z);

  double fNucleonTotalXsc, fNucleonElasticXsc, fNucleonInelasticXsc;
  double fNucleusTotalXsc, fNucleusElasticXsc, fNucleusInelasticXsc, fNucleusDiffractionXsc, fNucleusProductionXsc;
  
  static const double fNeutronBarCorrectionTot[93];
  static const double fNeutronBarCorrectionIn[93];
  
  static const double fProtonBarCorrectionTot[93];
  static const double fProtonBarCorrectionIn[93];

  static const double fPionPlusBarCorrectionTot[93];
  static const double fPionPlusBarCorrectionIn[93];

  static const double fPionMinusBarCorrectionTot[93];
  static const double fPionMinusBarCorrectionIn[93];

};

/////////////////////////////////////////////////////////////////////
//
// return correction at Tkin = 90*GeV GG -> Barashenkov tot xsc, when it 
// is available, else return 1.0


inline double HadronNucleusXsc::GetParticleBarCorTot(int particlePDG, int Z)
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

/////////////////////////////////////////////////////////////////////
//
// return correction at Tkin = 90*GeV GG -> Barashenkov in xsc, when it 
// is available, else return 1.0


inline double HadronNucleusXsc::GetParticleBarCorIn(int particlePDG, int Z)
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

#endif  // HadronNucleusXsc_H


// proton = 2212, neutron = 2112
// anti-proton = -2212
// pi+ = 211
// pi- = -211
// K+ = 321, K0L = 130
// K- = -321, K0S = 310
// sigma- = 3112, 
