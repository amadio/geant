#ifndef HepMCTruth_h
#define HepMCTruth_h

#include "MCTruthMgr.h"

#include "HepMC/Writer.h"
#include "HepMC/GenEvent.h"
#include "Geant/Typedefs.h"
#include "Geant/Config.h"
#include "GeantFwd.h"

#include <string>

namespace Geant {
  inline namespace GEANT_IMPL_NAMESPACE {
    class GeantTrack;
    struct MCEvent;
  }
}


namespace userapplication {

class HepMCTruth : public Geant::MCTruthMgr {
private:
  HepMC::Writer *output_file;

public:
  HepMCTruth();
  HepMCTruth(std::string &filename);
  ~HepMCTruth();

  virtual void InitMCTruthMgr();

  virtual bool CheckTrack(Geant::GeantTrack &gtrack, Geant::MCEvent* evt);
  
  virtual void CloseEvent(int evID);

  double fEMin; // minimum energy

private:
  HepMCTruth(const HepMCTruth &);            // no imp.
  HepMCTruth &operator=(const HepMCTruth &); // no imp.

};
 
}
#endif
