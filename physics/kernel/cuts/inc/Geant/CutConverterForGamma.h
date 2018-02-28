

#ifndef CUTCONVERTERFORGAMMA_H
#define CUTCONVERTERFORGAMMA_H

#include "Geant/CutConverter.h"

namespace geantphysics {
/**
 * @brief   Production threshold converter for gamma.
 * @class   CutConverterForGamma
 * @author  M Novak, A Ribon
 * @date    april 2016
 */
class CutConverterForGamma : public CutConverter {
public:
  CutConverterForGamma(int numebins = 301, double mincutenergy = 100.0 * geant::units::eV,
                       double maxcutenergy = 10.0 * geant::units::GeV);
  virtual ~CutConverterForGamma();

  virtual void Initialise();

protected:
  virtual void BuildLengthVector(const Material *mat);
  virtual double ComputeELossOrAbsXsecPerAtom(double zet, double ekin);

private:
  // some Z dependent cached variables for the approximated absorption cross section computation
  double fZ;
  double fS200keV;
  double fTmin;
  double fSmin;
  double fCmin;
  double fTlow;
  double fSlow;
  double fS1keV;
  double fClow;
  double fChigh;
};

} // namespace geantphysics

#endif // CUTCONVERTERFORGAMMA_H
