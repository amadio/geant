
#ifndef MyGVPhysicsList_h
#define MyGVPhysicsList_h 1

#include "G4VUserPhysicsList.hh"

class MyGVPhysicsList : public G4VUserPhysicsList {
public:
  MyGVPhysicsList();
  ~MyGVPhysicsList();

  virtual void ConstructParticle();
  virtual void ConstructProcess();
  //  virtual void SetCuts();

private:
  void BuildEMPhysics();
};

#endif
