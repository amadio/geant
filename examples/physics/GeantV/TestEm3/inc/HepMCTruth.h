#ifndef HepMCTruth_h
#define HepMCTruth_h

#include "Geant/MCTruthMgr.h"

#include "HepMC/Writer.h"
#include "HepMC/GenEvent.h"
#include "Geant/Typedefs.h"
#include "Geant/Config.h"
#include "Geant/Fwd.h"

#include <string>

namespace geant {
  inline namespace GEANT_IMPL_NAMESPACE {
    class Track;
    struct MCEvent;
  }
}


namespace userapplication {

class HepMCTruth : public geant::MCTruthMgr {
private:
  HepMC::Writer *output_file;

public:
  HepMCTruth();
  HepMCTruth(std::string &filename);
  ~HepMCTruth();

  virtual void InitMCTruthMgr();

  virtual bool CheckTrack(geant::Track &gtrack, geant::MCEvent* evt);
  
  virtual void CloseEvent(int evID);

  double fEMin; // minimum energy

private:
  HepMCTruth(const HepMCTruth &);            // no imp.
  HepMCTruth &operator=(const HepMCTruth &); // no imp.

};
 
}
#endif
