
#ifndef MyRun_h
#define MyRun_h 1

#include "G4Run.hh"
#include "globals.hh"
#include "MyRunDataPerPrimary.hh"

#include <vector>

class G4Track;
class MyEventDataPerPrimary;

class MyRun : public G4Run {
public:
  MyRun();
  virtual ~MyRun();

  virtual void Merge(const G4Run*);
  void FillPerEventNumPrimaries(G4int primtypeindx) { fRunDataPerPrimary[primtypeindx].fNumPrimaries += 1.;}
  void FillPerEvent(const MyEventDataPerPrimary&, G4int primtypeindx);
  void EndOfRun();

  const MyRunDataPerPrimary&  GetRunDataPerPrimary(G4int primtypeindx) const { return fRunDataPerPrimary[primtypeindx]; }


private:
  G4int   fNumPrimaryTypes;
  std::vector<MyRunDataPerPrimary>  fRunDataPerPrimary;

};

#endif
