
#ifndef ELASTICSCATTERINGPROCESS_H
#define ELASTICSCATTERINGPROCESS_H

#include "HadronicProcess.h"

#include "Proton.h"
#include "Neutron.h"
#include "PionPlus.h"
#include "PionMinus.h"
#include "PionZero.h"
#include "KaonPlus.h"
#include "KaonMinus.h"
#include "KaonZero.h"
#include "KaonShort.h"
#include "KaonLong.h"

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
