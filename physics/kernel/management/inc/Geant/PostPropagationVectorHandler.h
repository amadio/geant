
#ifndef POSTPROPAGATIONVECTORHANDLER_H
#define POSTPROPAGATIONVECTORHANDLER_H

// from geantV
#include "Geant/Handler.h"
// from geantV
namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {
class Propagator;
class Track;
class TaskData;
class Basket;
}
}

namespace geantphysics {

class MSCModel;

class PostPropagationVectorHandler : public geant::Handler {
public:
  /** @brief Default constructor */
  PostPropagationVectorHandler() : geant::Handler() {}

  /**
   * @brief Default constructor
   * @param propagator Propagator working with this handler
   */
  PostPropagationVectorHandler(int threshold, geant::Propagator *propagator, int modelIdx);

  /** @brief dtr */
  virtual ~PostPropagationVectorHandler();

  /** @brief Scalar DoIt interface */
  virtual void DoIt(geant::Track *track, geant::Basket &output, geant::TaskData *td);

  /** @brief Vector DoIt interface. Base class implements it as a loop. */
  virtual void DoIt(geant::Basket &input, geant::Basket &output, geant::TaskData *td);

private:
  PostPropagationVectorHandler(const PostPropagationVectorHandler &) = delete;
  PostPropagationVectorHandler &operator=(const PostPropagationVectorHandler &) = delete;

  MSCModel *fModel;
};

} // namespace geantphysics

#endif // POSTPROPAGATIONHANDLER_H
