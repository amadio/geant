#ifndef NeutronNudyInelasticXsec_H
#define NeutronNudyInelasticXsec_H

#include "Geant/HadronicCrossSection.h"
#include "Geant/TNudyEndfRecoPoint.h"
#include <string>

namespace geantphysics {
class NeutronNudyInelasticXsec : public HadronicCrossSection {
public:
  NeutronNudyInelasticXsec();
  virtual ~NeutronNudyInelasticXsec();

  double GetIsotopeCrossSection(const int particleCode, const double energyKin, const double mass, 
                                const int Z, const int A);

private:
  const char *fRENDF;
  std::string filename;
  NudyPhysics::TNudyEndfRecoPoint recopoint;
};

}
#endif // NeutronNudyInelasticXsec_H
