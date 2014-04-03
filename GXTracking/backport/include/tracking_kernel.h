#ifndef TRACKING_GPU_H
#define TRACKING_GPU_H

#include "G4VPhysicalVolume.hh"
#include "GPFieldMap.h"
#include "G4Track.hh"

#include "G4PhysicsTable.hh"
#include <cmath>

class G4MagneticField;
class G4Navigator;

extern void tracking_cpu(G4VPhysicalVolume *world, G4Navigator &aNavigator,
			 G4MagneticField *magMap, 
			 G4Track **track, 
			 G4PhysicsTable* eBrem_table, 
			 G4PhysicsTable* eIoni_table, 
			 G4PhysicsTable* msc_table, 
			 size_t nTrackSize); 

extern void tracking_electron_cpu(G4VPhysicalVolume *world, G4Navigator &aNavigator,
				  G4MagneticField *magMap, 
				  G4Track **track, 
				  G4PhysicsTable* eBrem_table, 
				  G4PhysicsTable* eIoni_table, 
				  G4PhysicsTable* msc_table, 
				  size_t nTrackSize); 

extern void tracking_photon_cpu(G4VPhysicalVolume *world, G4Navigator &aNavigator,
				G4MagneticField *magMap, 
				G4Track **track, 
				G4PhysicsTable* eBrem_table, 
				G4PhysicsTable* eIoni_table, 
				G4PhysicsTable* msc_table, 
				size_t nTrackSize); 

#endif
