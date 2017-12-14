#ifndef HepMCTruth_h
#define HepMCTruth_h

#include "MCTruthMgr.h"

#include "HepMC/Writer.h"
#include "HepMC/GenEvent.h"
#include "GeantFwd.h"

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

#ifdef USE_ROOT
  ClassDef(HepMCTruth, 1)
#endif
};

#endif
