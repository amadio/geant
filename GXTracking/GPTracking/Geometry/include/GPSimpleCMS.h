// Renamed CMSModel.hpp created by Otto Seiskari, August 2010
// Parses a subset of geometry of the CMS detector from cms.txt
// ENABLE_VOXEL_NAVIGATION=1
// ENABLE_SLICED_TUBS=1
// ENABLE_G4POLYCONE=0

#ifndef GPSimpleCMS_HH
#define GPSimpleCMS_HH

#include "GPUserGeometry.h"
#include "GPBox.h"
#include "GPTubs.h"
#include "GPCons.h"

#include <fstream>
#include <map>
#include <list>

class GPSimpleCMS : public GPUserGeometry
{
public:
  GPSimpleCMS( int maxsz = 1024*1000000 ) : 
    GPUserGeometry ( maxsz ),                // ~1GB
    scalefactor ( 1E-3 )                     // Scale in meters!
  {}
	
  double getScale() const
  {
    return 28600*scalefactor;
  }

  /*
  CameraParameters getCamera() const
  {
    CameraParameters p;
    p.heading = 20;
    p.pitch = 60;
    p.dist = 0.22*50000*scalefactor;
    p.yfov = 80;
    p.target_z = 0.04*50000*scalefactor;
    return p;
  }
  */
  
private:

  const G4double scalefactor;

  GPVPhysicalVolume *worldVolume;
  GPLogicalVolume *logWorld;

  std::map<std::string, GPMaterial*> materials;
  std::map<std::string, GPVSolid*> solids;
  std::map<std::string, GPLogicalVolume*> logvols;
  std::map<GPVPhysicalVolume*, std::string> phystologref;
  std::list<GPVPhysicalVolume*> physvols;

  void create_impl()
  {
    const char *filename = "./cms.txt";
    std::ifstream stream(filename);
    
    worldVolume = reserveThing<GPVPhysicalVolume>();
    
    readMaterials( stream );
    readSolids( stream );
    readStructure( stream );

    //#ifdef VERBOSE
    std::cerr << "Imported " << materials.size() << " materials\n";
    std::cerr << "Imported " << solids.size() << " solids\n";
    std::cerr << "Imported " << logvols.size() << " logical and "
	      << physvols.size()+1 << " physical volumes\n";
    //#endif
		
    for ( std::list<GPVPhysicalVolume*>::const_iterator itr = physvols.begin();
	  itr != physvols.end(); ++itr ) {
      addPhysicalVolumePointers( *itr );
    }

    for ( std::map<std::string,GPLogicalVolume*>::const_iterator itr = 
	    logvols.begin(); itr != logvols.end(); ++itr ) {
      //#ifdef ENABLE_VOXEL_NAVIGATION
      /*
      const int VOX_THRESHOLD = GPVoxelHeader_GetMinVoxelVolumesLevel1();
      const int nd = GPLogicalVolume_GetNoDaughters(itr->second);
      if ( nd >= VOX_THRESHOLD )
      {
	createVoxels( itr->second );
      }
      */
      //#endif
      addLogicalVolumePointers( itr->second );
    }
  }
	
  template <class T> static T get( std::istream& in )
  {
    T obj;
    in >> obj;
    if (!in) throw std::runtime_error("Failed to read input");
    return obj;
  }
  
  static G4double degToRad( G4double d )
  {
    return d / 360.0 * twopi;
  }

  GPThreeVector readVec( std::istream& in )
  {
    const G4double x = get<G4double>(in)*scalefactor;
    const G4double y = get<G4double>(in)*scalefactor;
    const G4double z = get<G4double>(in)*scalefactor;
    return GPThreeVector_create(x,y,z);
  }
	
  GPRotationMatrix readRot( std::istream& in )
  {
    G4double elem[3][3];
    for (int i=0; i<3; ++i)
      for (int j=0; j<3; ++j)
	elem[i][j] = get<G4double>(in);
    
    return GPRotationMatrix_create(elem[0][0], elem[0][1], elem[0][2],
				   elem[1][0], elem[1][1], elem[1][2],
				   elem[2][0], elem[2][1], elem[2][2]);
  }
	
  void readMaterials( std::istream& in )
  {
    const int numberOfMaterials = get<int>( in );
    
    // Possibly override material settings
    
    const enum { OVERRIDE, INCLUDE, EXCLUDE, NONE }
    overrideMaterials = INCLUDE;
    
    const G4double default_property = 1.0;
    
    // The densities are given as g/cm^3, lengths are already scaled to meters
    // multiply by 1000 to get the densities in kg/m^3
    const G4double densityMultiplier = 1000.0;
    
    std::map<std::string, G4double> materialProperties;
    // Leave empty, no overrides
    
    GPMaterial *matarray = reserveNThings<GPMaterial>(numberOfMaterials);
    
    for ( int i=0; i<numberOfMaterials; ++i ) {
      const std::string materialName = get<std::string>(in);
      const G4double D = get<G4double>(in);
      const G4double atom = get<G4double>(in);
      // Take "density" value from "D" value of material
      G4double property = D + 0*atom;	
      
      if ( materialProperties.count( materialName ) ) {
	switch ( overrideMaterials ) {
	case EXCLUDE: property = default_property; break;
	case OVERRIDE: property = materialProperties[materialName]; break;
	default: break;
	}
      }
      else {
	switch ( overrideMaterials ) {
	case INCLUDE: property = default_property; break;
	default: break;
	}
      }
      
      matarray[i].fDensity = property*densityMultiplier;
      
      if (materials.count(materialName))
	throw std::runtime_error("Duplicate material key "+materialName);
      
      materials[materialName] = matarray+i;
    }
  }
  
  void readSolids( std::istream& in )
  {
    const int numberOfSolids = get<int>( in );
    
    for ( int i=0; i<numberOfSolids; ++i ) {
      const std::string type = get<std::string>(in);
      const std::string name = get<std::string>(in);
      
      // std::cerr << "Importing " << type << " " << name << "\n";
      
      GPVSolid *s = NULL;
      
      if ( type == "Box" ) {
	const GPThreeVector xyz = readVec(in);
	GPBox *box = reserveThing<GPBox>();
	GPBox_Constructor( box, xyz.x, xyz.y, xyz.z );
	s = (GPVSolid*)box;
      }
      else if ( type == "Tube" ) {
	const G4double rmin = get<G4double>(in)*scalefactor;
	const G4double rmax = get<G4double>(in)*scalefactor;
	const G4double z = get<G4double>(in)*scalefactor;
	const G4double sphi = degToRad(get<G4double>(in));
	const G4double dphi = degToRad(get<G4double>(in));
	
	/*
        #ifndef ENABLE_SLICED_TUBS
	if (std::fabs(dphi-twopi) > 0.1) {
	  G4Cons *tubs = reserveThing<G4Cons>();
	  s = (GPVSolid*)tubs;
	  G4Cons_ctor( tubs, rmin, rmax, rmin, rmax, z, sphi, dphi );
	}
	else {
        #endif
	*/
	GPTubs *tubs = reserveThing<GPTubs>();
	s = (GPVSolid*)tubs;
	GPTubs_Constructor( tubs, rmin, rmax, z, sphi, dphi );
	//#ifndef ENABLE_SLICED_TUBS
	//	}
	//#endif
      }
      else if ( type == "Cone" ) {
	GPCons *cons = reserveThing<GPCons>();
	s = (GPVSolid*)cons;
	    
	const G4double rmin1 = get<G4double>(in)*scalefactor;
	const G4double rmax1 = get<G4double>(in)*scalefactor;
	const G4double rmin2 = get<G4double>(in)*scalefactor;
	const G4double rmax2 = get<G4double>(in)*scalefactor;
	const G4double z = get<G4double>(in)*scalefactor;
	const G4double sphi = degToRad(get<G4double>(in));
	const G4double dphi = degToRad(get<G4double>(in));
	
	GPCons_Constructor( cons, rmin1, rmax1, rmin2, rmax2, z, sphi, dphi );
      }
      /*
      #ifdef ENABLE_G4POLYCONE
      else if ( type == "Polycone" ) {
	const G4double sphi = degToRad(get<G4double>(in));
	const G4double dphi = degToRad(get<G4double>(in));
	const int nplanes = get<int>( in );
	    
	G4double *zarr = reserveNThings<G4double>(nplanes);
	G4double *rminarr = reserveNThings<G4double>(nplanes);
	G4double *rmaxarr = reserveNThings<G4double>(nplanes);
	    
	for ( int j=0; j<nplanes; ++j ) {
	  if ( get<std::string>(in) != "ZPlane" )
	    throw std::runtime_error( "Invalid input" );
		
	  const G4double r2 = get<G4double>(in)*scalefactor;
	  const G4double r1 = get<G4double>(in)*scalefactor;
	  const G4double z = get<G4double>(in)*scalefactor;
	  
	  zarr[j] = z;
	  rminarr[j] = r1;
	  rmaxarr[j] = r2;
	}
	    
	G4PolyCone *cons = reserveThing<G4PolyCone>();
	s = (GPVSolid*)cons;
	
	G4PolyCone_ctor( cons, nplanes, zarr, rminarr, rmaxarr, sphi, dphi );
	addPolyConePointers( cons );
      }
      #endif
      */
      else
	throw std::runtime_error("Unsupported object type: "+type);
      assert( s != NULL );
      
      if (solids.count(name))
	throw std::runtime_error("Duplicate solid key "+name);
      
      solids[ name ] = s;
    }
  }
  
  
  void readStructure( std::istream &in )
  {
    const int numberOfLogicalVolumes = get<int>(in);
    
    for ( int i=0; i<numberOfLogicalVolumes; ++i ) {
      const std::string volName = get<std::string>(in);
      const std::string materialKey = get<std::string>(in);
      const std::string solidKey = get<std::string>(in);
      
      GPLogicalVolume *logVol = reserveThing<GPLogicalVolume>();		
      GPLogicalVolume_Constructor( logVol,
				   solids.find(solidKey)->second,
				   materials.find(materialKey)->second );
	
      if (logvols.count(volName))
	throw std::runtime_error("Duplicate volume key "+volName);
	
      logvols[volName] = logVol;
	
      if ( i == 0 ) {
	logWorld = logVol;
	const GPRotationMatrix idRot = GPRotationMatrix_create(1,0,0,
							       0,1,0,
							       0,0,1);
	const GPThreeVector zeroTrans = GPThreeVector_create(0,0,0);
	GPVPhysicalVolume_Constructor( worldVolume, idRot, zeroTrans, logVol );
	addPhysicalVolumePointers( worldVolume );
      }
	
      const int numDaughters = get<int>(in);
      if ( numDaughters > 0 ) {
	GPVPhysicalVolume *pdarr =
	  reserveNThings<GPVPhysicalVolume>(numDaughters);
	
	for ( int j=0; j<numDaughters; ++j ) {
	  if ( get<std::string>(in) != "PhysicalVolume" )
	    throw std::runtime_error("Invalid input");
	  // ignore 'PhysicalVolume'
	  
	  GPVPhysicalVolume *physvol = pdarr+j;
	  phystologref[physvol] = get<std::string>(in);
		
	  const GPRotationMatrix rot = readRot(in);
	  const GPThreeVector trans = readVec(in);
		
	  GPVPhysicalVolume_Constructor( physvol, rot, trans, NULL );
	  GPVPhysicalVolume_SetMotherLogical( physvol, logVol );
	  physvols.push_back(physvol);
	  addLogicalVolumeDaughter( logVol, physvol );
	}
      }
    }
    
    for ( std::list<GPVPhysicalVolume*>::const_iterator itr = physvols.begin();
	  itr != physvols.end(); ++itr ) {
      GPVPhysicalVolume_SetLogicalVolume( *itr,
      logvols.find(phystologref.find(*itr)->second)->second );
    }
  }
};

#endif
