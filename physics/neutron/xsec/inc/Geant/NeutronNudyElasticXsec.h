#ifndef NeutronNudyElasticXsec_H
#define NeutronNudyElasticXsec_H

#include "Geant/HadronicCrossSection.h"
#include "Geant/TNudyEndfRecoPoint.h"
#include <string>

namespace geantphysics {
class NeutronNudyElasticXsec : public HadronicCrossSection {
public:
  NeutronNudyElasticXsec();
  virtual ~NeutronNudyElasticXsec();

  double GetIsotopeCrossSection(const int particleCode, const double energyKin, const double mass, const int Z,
                                const int A);

private:
  const char *fRENDF;
  std::string filename;
  NudyPhysics::TNudyEndfRecoPoint recopoint;
};

} // namespace geantphysics
#endif // NeutronNudyElasticXsec_H
