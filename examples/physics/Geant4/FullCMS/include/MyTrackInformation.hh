
#ifndef MyTrackInformation_h
#define MyTrackInformation_h 1

#include "G4VUserTrackInformation.hh"
#include "globals.hh"

class MyTrackInformation : public G4VUserTrackInformation {
public:
  MyTrackInformation(G4int id);
  virtual ~MyTrackInformation();

  G4int GetPrimaryTrackID() const { return fPrimaryTrackID; }
  void SetPrimaryTrackID(G4int id) { fPrimaryTrackID = id; }

private:
  G4int fPrimaryTrackID;
};

#endif
