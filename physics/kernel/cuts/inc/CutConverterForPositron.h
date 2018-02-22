#ifndef CUTCONVERTERFORPOSITRON_H
#define CUTCONVERTERFORPOSITRON_H

#include "CutConverter.h"

namespace geantphysics {
/**
 * @brief   Production threshold converter for positron.
 * @class   CutConverterForPositron
 * @author  M Novak, A Ribon
 * @date    april 2016
 */
class CutConverterForPositron : public CutConverter {
public:
  CutConverterForPositron(int numebins = 301, double mincutenergy = 100.0*geant::units::eV,
                          double maxcutenergy = 10.0*geant::units::GeV);
  virtual ~CutConverterForPositron();

  virtual void   Initialise();

protected:
  virtual double ComputeELossOrAbsXsecPerAtom(double zet, double ekin);

};

} // namespace geantphysics

#endif // CUTCONVERTERFORPOSITRON_H
