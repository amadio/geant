#ifndef GPSimpleEcal_HH
#define GPSimpleEcal_HH

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4PVPlacement.hh"

#include "GPUserGeometry.h"

class GPSimpleEcal : public GPUserGeometry
{
public:
 GPSimpleEcal( int nphi, int nz, double density)
    : GPUserGeometry( 1024 * nphi * nz ), fWorldLog(0), fWorldPhy(0) 
  {
    fNphi = nphi;
    fNz = nz;
    fDensity = density;
  }

  double getScale() const { return 1.0; }

private:
  int fNphi;
  int fNz;
  double fDensity;
  G4LogicalVolume   *fWorldLog;
  G4VPhysicalVolume *fWorldPhy;

private:

  void create_impl()
  {
    // World volume - CMS dimension

    const G4double world_x =  9000.;
    const G4double world_y =  9000.;
    const G4double world_z = 16000;

    G4RotationMatrix *idRot = new G4RotationMatrix();/* (1,0,0, */
                                                   /* 0,1,0, */
                                                   /* 0,0,1); */
    
    G4ThreeVector zeroTrans( 0,0,0 );

    G4VPhysicalVolume *world_phy = 0; // reserveThing<GPVPhysicalVolume>();
    G4LogicalVolume   *world_log = 0; // reserveThing<GPLogicalVolume>();
    G4Box *world =  0; // reserveThing<GPBox>();
    G4Material *world_mat = 0; // reserveThing<GPMaterial>();
    //world_mat->fDensity = 0.0;

    world = new G4Box("world", world_x, world_y, world_z);
    const G4double universe_mean_density = 1.e-25*g/cm3;
    world_mat = new G4Material("Galactic", 1,1.01*g/mole , universe_mean_density); // world_mat->fDensity,1,2);
    world_log = new G4LogicalVolume(world, world_mat, "world_log");
    world_phy = new G4PVPlacement(idRot, zeroTrans, "simpleEcal", world_log, 0, false, 0);

    // Construct Calorimeter

    int crystal_nphi  = fNphi;
    int crystal_nz    = fNz;
    G4double ecal_density = fDensity;

    const int crystal_n     = crystal_nphi*crystal_nz;
    const G4double ecal_zmin  = -3000.;
    const G4double ecal_zmax  =  3000.;

    const G4double ecal_rmin  =  10.;
    const G4double ecal_rmax  =  5000.;
    const G4double ecal_dz    =  0.5*(ecal_zmax-ecal_zmin)/crystal_nz;
    const G4double ecal_sphi  =     0.;
    const G4double ecal_dphi  =  2.0*M_PI/crystal_nphi;

    G4VPhysicalVolume **ecal_phy = new G4VPhysicalVolume*[crystal_n]; // reserveNThings<GPVPhysicalVolume>(crystal_n);
    G4LogicalVolume   **ecal_log = new G4LogicalVolume*[crystal_n]; // reserveNThings<GPLogicalVolume>(crystal_n);

    int iptr = 0;

    for ( int j = 0; j < crystal_nz ; ++j ) {
      for ( int i = 0; i < crystal_nphi ; ++i ) {

        iptr = i+j*crystal_nphi;
        // G4VSolid *solid = NULL;

	G4Tubs *ecal = 0; // reserveThing<GPTubs>();
	G4Material *ecal_mat = 0; //reserveThing<GPMaterial>();
	
   G4Element* elPb = new G4Element( "Lead", "Pb", 82., 207.19*g/mole );
   G4Element* elW = new G4Element( "Tungstenm", "W",74., 183.85*g/mole);
   G4Element* elO = new G4Element( "Oxygen", "O2", 8., 16.*g/mole );

         
   ecal = new G4Tubs("tubs",
                          ecal_rmin, ecal_rmax,  
                          ecal_dz, ecal_sphi+i*ecal_dphi,ecal_dphi);
	ecal_mat = // new G4Material("FakeStuff",90,120,ecal_density); // ecal_mat,ecal_density,90,120);
     new G4Material("PbWO4", 8.28*g/cm3, 3);
   ecal_mat ->AddElement( elPb, 1 );
   ecal_mat ->AddElement( elW,  1 );
   ecal_mat ->AddElement( elO,  4 );

	G4ThreeVector ecal_trans(0,0,ecal_zmin+(2.0*j+1.0)*ecal_dz);

	ecal_log[iptr] = new G4LogicalVolume(ecal, ecal_mat, "crystal");
	ecal_phy[iptr] = new G4PVPlacement(idRot, ecal_trans, "crystal", ecal_log[iptr], 
                                           world_phy, //0, 
                                           false, 0);
	
	//Set mother
	ecal_phy[iptr]->SetMotherLogical(world_log);

	addLogicalVolumePointers( ecal_log[iptr]);
	addPhysicalVolumePointers( ecal_phy[iptr]);

      }      
    }

    //add daughter volume 
    for ( int j=0; j < crystal_nz ; ++j ) {
      for ( int i=0; i < crystal_nphi ; ++i ) {
        iptr = i+j*crystal_nphi;
	addLogicalVolumeDaughter( world_log, ecal_phy[iptr]);
      }    
    }

    // Register world volume pointers for relocation
    addLogicalVolumePointers( world_log );
    addPhysicalVolumePointers( world_phy );
    fWorldLog = world_log;
    fWorldPhy = world_phy;

  }

  // Get the world volume.
  G4VPhysicalVolume *getWorldVolume() const
  {
     return fWorldPhy;
  }

  // Get the world volume.
  G4LogicalVolume *getWorldLogicalVolume() const
  {
     return fWorldLog;
  }

};

#endif
