//
// Differential cross section class for Compton scattering using Klein-Nishina formulation
// 
// Independent output variable:  Energy of outgoing photon
// 
#include "GUTypeDef.h"

class GUXSectionKleinNishina
{
  public: 
    FQUALIFIER GUXSectionKleinNishina(){}
    FQUALIFIER ~GUXSectionKleinNishina(){}

    FQUALIFIER double CalculateDiffCrossSection( int Zelement, double Ein, double outEphoton ) const;


private:
	// no state since this class is calculating the diff cross section from a formula


};

// function implementing the cross section for KleinNishina
// TODO: need to get electron properties from somewhere
const double electron_mass_c2 = 0.51;

FQUALIFIER
double CalculateDiffCrossSection( int Zelement, double energy0, 
				  double energy1 )
{
  // based on Geant4 : G4KleinNishinaCompton
  // input  : energy0 (incomming photon energy)
  //          energy1 (scattered photon energy)
  // output : dsigma  (differential cross section) 

  double E0_m = energy0/electron_mass_c2 ;
  double epsilon = energy1/energy0;

  double onecost = (1.- epsilon)/(epsilon*E0_m);
  double sint2   = onecost*(2.-onecost);
  double greject = 1. - epsilon*sint2/(1.+ epsilon*epsilon);
  double dsigma = (epsilon + 1./epsilon)*greject;

  return dsigma;
}
