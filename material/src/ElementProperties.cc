
#include "ElementProperties.h"

#include "PhysicalConstants.h"
#include "Element.h"

#include <cmath>

namespace geantphysics {

ElementProperties::ElementProperties(Element *elem) : fElement(elem) {
  InitialiseMembers();
}

ElementProperties::~ElementProperties() {}

void ElementProperties::InitialiseMembers() {
  double z = fElement->GetZ();
  fZ13    = std::pow(z,1./3.);
  fZ23    = std::pow(z,2./3.);
  fLogZ   = std::log(z);
  fLogZ13 = fLogZ/3.;
  fLogZ23 = 2.*fLogZ/3.;
  ComputeCoulombCorrection(z);
}

/**
 * @internal
 *
 *  The Coulomb correction \f$f\f$ is computed as \cite davies1954theory [Eqs.(36-38)]
 *  \f[
 *   f(\nu) = \nu^2 \sum_{n=1}^{\infty} \frac{1}{n(n^2+\nu^2)} = \nu^2 \left[  1/(1+\nu^2)  + 0.20206 - 0.0369\nu^2
 *            + 0.0083\nu^4 - 0.002\nu^6 \right]
 *  \f]
 * where \f$\nu=\alpha Z\f$ with \f$ \alpha \f$ being the fine structure constant.
 *
 * @endinternal
 */
void ElementProperties::ComputeCoulombCorrection(double z) {
  // Coulomb correction from Davis,Bethe,Maximom PRL 1954 Eqs.(36-38)
  double mu  = z*geant::kFineStructConst;
  double mu2 = mu*mu;
  double mu4 = mu2*mu2;
  double mu6 = mu2*mu4;
  fCoulombCorrection = mu2*(1.0/(1.0+mu2) + 0.20206 - 0.0369*mu2 + 0.0083*mu4 - 0.002*mu6);
}

} // namespace geantphysics
