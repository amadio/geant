#ifndef NeutronNudyFissionXsec_H
#define NeutronNudyFissionXsec_H

#include "Geant/HadronicCrossSection.h"
#include "Geant/TNudyEndfRecoPoint.h"
#include <string>

namespace geantphysics {
class NeutronNudyFissionXsec : public HadronicCrossSection {
public:
  NeutronNudyFissionXsec();
  virtual ~NeutronNudyFissionXsec();

  double GetIsotopeCrossSection(const int particleCode, const double energyKin, const double mass, 
                                const int Z, const int A);

private:
  const char *fRENDF;
  std::string filename;
  NudyPhysics::TNudyEndfRecoPoint recopoint;
};

}
#endif // NeutronNudyFissionXsec_H
