#ifndef NeutronNudyCaptureXsec_H
#define NeutronNudyCaptureXsec_H

#include "Geant/HadronicCrossSection.h"
#include "Geant/TNudyEndfRecoPoint.h"
#include <string>

namespace geantphysics {
class NeutronNudyCaptureXsec : public HadronicCrossSection {
public:
  NeutronNudyCaptureXsec();
  virtual ~NeutronNudyCaptureXsec();

  double GetIsotopeCrossSection(const int particleCode, const double energyKin, const double mass, const int Z,
                                const int A);

private:
  const char *fRENDF;
  std::string filename;
  NudyPhysics::TNudyEndfRecoPoint recopoint;
};

} // namespace geantphysics
#endif // NeutronNudyCaptureXsec_H
