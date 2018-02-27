#ifndef Parameterizations_H
#define Parameterizations_H

namespace geantphysics {

/// Mandelstam
double CalcMandelstamS(const double mp, const double mt, const double Plab);
/// Nucleus Radius
double GetNucleusRadius(int At);
/// Coulomb barrier
double GetCoulombBarrier(int particlePDG, double proj_mass, double energyKin, int targetPDG, double target_mass);
//
// Returns hadron-nucleon Xsc according to PDG parametrisation (2005):
// http://pdg.lbl.gov/2006/reviews/hadronicrpp.pdf
double GetHadronNucleonXscPDG(int particlePDG, double mass, double energyKin, int targetPDG);
//
// Returns hadron-nucleon cross-section based on N. Starkov parametrisation of
// data from mainly http://wwwppds.ihep.su:8001/c5-6A.html database
double GetHadronNucleonTotalXscNS(int particlePDG, double mass, double energyKin, int targetPDG);
//
// Returns kaon-nucleon cross-section based on smoothed NS for GG model
double GetKaonNucleonTotalXscGG(int particlePDG, double mass, double energyKin, int targetPDG);
//
// Returns hadron-nucleon inelastic cross-section based on N. Starkov parametrisation of
// data from mainly http://wwwppds.ihep.su:8001/c5-6A.html database
double GetHadronNucleonInelasticXscNS(int particlePDG, double mass, double energyKin, int targetPDG);
//
// Returns kaon-nucleon inelastic cross-section based on smoothed NS for GG model
double GetKaonNucleonInelasticXscGG(int particlePDG, double mass, double energyKin, int targetPDG);

}      // namespace geantphysics

#endif // Parameterizations_H
