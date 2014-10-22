//
// Differential cross section class for Compton scattering using Klein-Nishina formulation
// 
// Independent output variable:  Energy of outgoing photon
// 
class GUXSectionKleinNishina
{
  public: 
    FQUALIFIER GUXSectionKleinNishina();
    FQUALIFIER ~GUXSectionKleinNishina();

    FQUALIFIER double CalculateDiffCrossSection( int Zelement, double Ein, double outEphoton );
};