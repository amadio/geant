#ifndef GXTRACK_H
#define GXTRACK_H 

struct GXTrack
{
  G4int status;
  G4int id;
  G4int proc;
  G4double x; 
  G4double y;
  G4double z;
  G4double px;
  G4double py;
  G4double pz;
  G4double E;
  G4double q;
  G4double s;
} ;

#endif
