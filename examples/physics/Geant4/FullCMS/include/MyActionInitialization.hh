
#ifndef MyActionInitialization_h
#define MyActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class MyActionInitialization : public G4VUserActionInitialization {
public:
  MyActionInitialization(bool isperformance = false);
  virtual ~MyActionInitialization();

  virtual void BuildForMaster() const;
  virtual void Build() const;

  void SetPerformanceModeFlag(bool val) { fIsPerformance = val; }

private:
  bool fIsPerformance;
};

#endif
