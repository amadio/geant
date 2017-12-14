#ifndef GlauberGribovElasticXsc_H
#define GlauberGribovElasticXsc_H

#include "HadronicCrossSection.h"
#include "GlauberGribovTotalXsc.h"
#include "GlauberGribovInelasticXsc.h"

namespace geantphysics {

  class GlauberGribovElasticXsc : public HadronicCrossSection
  {
  public:

    GlauberGribovElasticXsc ();
    virtual ~GlauberGribovElasticXsc ();
   
    double GetIsotopeCrossSection(const int particleCode, const double energyKin, const double mass,
				  const int Z, const int N);
    
  private:

    GlauberGribovTotalXsc* GGTotalXsc;
    GlauberGribovInelasticXsc* GGInelasticXsc;
  };

}       // namespace geantphysics

#endif  // GlauberGribovElasticXsc_H

