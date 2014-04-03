// Renamed CMSModel.hpp created by Otto Seiskari, August 2010
// Parses a subset of geometry of the CMS detector from cms.txt
// ENABLE_VOXEL_NAVIGATION=1
// ENABLE_SLICED_TUBS=1
// ENABLE_G4POLYCONE=0

#ifndef GXSimpleCMS_HH
#define GXSimpleCMS_HH

#include "GXUserGeometry.h"
#include <fstream>
#include <map>
#include <list>

class GXSimpleCMS : public GXUserGeometry
{
public:
  GXSimpleCMS( int maxsz = 1024*1000000 ) : 
    GXUserGeometry ( maxsz ),                // ~1GB
    scalefactor ( 1E-3 )                     // Scale in meters!
  {}
	
  double getScale() const
  {
    return 28600*scalefactor;
  }

private:

  const G4double scalefactor;

  GPVPhysicalVolume *worldVolume;
  GPLogicalVolume *logWorld;

  std::map<std::string, GPMaterial*> materials;
  std::map<std::string, GPVSolid*> solids;
  std::map<std::string, GPLogicalVolume*> logvols;
  std::map<GPVPhysicalVolume*, std::string> phystologref;
  std::list<GPVPhysicalVolume*> physvols;

  void create_impl();
	
  template <class T> static T get( std::istream& in );
  
  static G4double degToRad( G4double d );

  GPThreeVector readVec( std::istream& in );
	
  GPRotationMatrix readRot( std::istream& in );
	
  void readMaterials( std::istream& in );
  
  void readSolids( std::istream& in );
  
  void readStructure( std::istream &in );

};

#endif
