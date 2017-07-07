
#ifndef ELASTICSCATTERINGPROCESS_H
#define ELASTICSCATTERINGPROCESS_H

#include "HadronicProcess.h"

#include "Proton.h"
#include "Neutron.h"

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
