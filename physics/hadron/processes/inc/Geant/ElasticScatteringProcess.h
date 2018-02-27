
#ifndef ELASTICSCATTERINGPROCESS_H
#define ELASTICSCATTERINGPROCESS_H

#include "Geant/HadronicProcess.h"

#include "Geant/Proton.h"
#include "Geant/Neutron.h"
#include "Geant/PionPlus.h"
#include "Geant/PionMinus.h"
#include "Geant/PionZero.h"
#include "Geant/KaonPlus.h"
#include "Geant/KaonMinus.h"
#include "Geant/KaonZero.h"
#include "Geant/KaonShort.h"
#include "Geant/KaonLong.h"

#include <string>

namespace geantphysics {

/**
 * @brief   Hadron elastic scattering process.
 * @class   ElasticScatteringProcess
 * @author  W Pokorski
 * @date    June 2017
 */


class ElasticScatteringProcess : public HadronicProcess {
public:
  // CTR
  ElasticScatteringProcess(const std::string &name = "hElastic");
};

} // namespace geantphysics

#endif
