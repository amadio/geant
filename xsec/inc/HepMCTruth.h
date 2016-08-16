#ifndef HepMCTruth_h
#define HepMCTruth_h

#include "MCTruthMgr.h"

#if __cplusplus >= 201103L
#include "HepMC/Writer.h"
#include "HepMC/GenEvent.h"
#endif

class HepMCTruth : public MCTruthMgr {
private:
#if __cplusplus >= 201103L
  HepMC::Writer *output_file;
#endif

public:
  HepMCTruth();
  HepMCTruth(std::string &filename);
  ~HepMCTruth();

  virtual void InitMCTruthMgr();

  virtual bool CheckTrack(Geant::GeantTrack &gtrack, MCEvent* evt);
  
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
