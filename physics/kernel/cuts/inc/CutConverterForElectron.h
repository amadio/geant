
#ifndef CUTCONVERTERFORELECTRON_H
#define CUTCONVERTERFORELECTRON_H

#include "CutConverter.h"

namespace geantphysics {
/**
 * @brief   Production threshold converter for electron.
 * @class   CutConverterForElectron
 * @author  M Novak, A Ribon
 * @date    april 2016
 */
class CutConverterForElectron : public CutConverter {
public:
  CutConverterForElectron(int numebins = 301, double mincutenergy = 100.0*geant::eV,
                          double maxcutenergy = 10.0*geant::GeV);
  virtual ~CutConverterForElectron();

  virtual void   Initialise();

protected:
  virtual double ComputeELossOrAbsXsecPerAtom(double zet, double ekin);

};

} // namespace geantphysics

#endif // CUTCONVERTERFORELECTRON_H
