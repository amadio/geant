#ifndef GPSimpleEcal_HH
#define GPSimpleEcal_HH

#include "GPBox.h"
#include "GPTubs.h"
#include "GPThreeVector.h"
#include "GPRotationMatrix.h"

#include "GPUserGeometry.h"
#include "GPGeomManager.h"

class GPSimpleEcal : public GPUserGeometry
{
public:
   GPSimpleEcal( int nphi, int nz, double density)
   : GPUserGeometry( 1024 * nphi * nz )
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
   
private:
   
   void create_impl()
   {
      // World volume - CMS dimension
      
      const G4double world_x =  9000.;
      const G4double world_y =  9000.;
      const G4double world_z = 16000;
      
      GPRotationMatrix idRot = GPRotationMatrix_create(1,0,0,
                                                       0,1,0,
                                                       0,0,1);
      
      GPThreeVector zeroTrans = GPThreeVector_create( 0,0,0 );
      
      GPGeomManager *geomManager = reserveThing<GPGeomManager>();
      
      GPVPhysicalVolume *world_phy = reserveThing<GPVPhysicalVolume>();
      
      geomManager->fWorldPhysical = world_phy;
      addPointer(geomManager->fWorldPhysical);
      
      GPLogicalVolume   *world_log = reserveThing<GPLogicalVolume>();
      GPBox *world = reserveThing<GPBox>();
      GPMaterial *world_mat = reserveThing<GPMaterial>();
      world_mat->fDensity = 0.0;
      
      GPBox_Constructor(world, world_x, world_y, world_z);
      GPMaterial_Constructor(world_mat,world_mat->fDensity,1,2);
      GPLogicalVolume_Constructor(world_log, (GPVSolid*)world, world_mat);
      GPVPhysicalVolume_Constructor(world_phy, idRot, zeroTrans, world_log);
      
      // Construct Calorimeter
      
      int crystal_nphi  = fNphi;
      int crystal_nz    = fNz;
      // G4double ecal_density = fDensity;
      
      const int crystal_n     = crystal_nphi*crystal_nz;
      const G4double ecal_zmin  = -3000.;
      const G4double ecal_zmax  =  3000.;
      
      const G4double ecal_rmin  =  10.;
      const G4double ecal_rmax  =  5000.;
      const G4double ecal_dz    =  0.5*(ecal_zmax-ecal_zmin)/crystal_nz;
      const G4double ecal_sphi  =     0.;
      //    const G4double ecal_dphi  =  2.0*M_PI/crystal_nphi;
      const G4double ecal_dphi  =  2.0*pi/crystal_nphi;
      
      GPVPhysicalVolume *ecal_phy = reserveNThings<GPVPhysicalVolume>(crystal_n);
      GPLogicalVolume   *ecal_log = reserveNThings<GPLogicalVolume>(crystal_n);
      
      int iptr = 0;

      //Build Material: PbWO4 a.k.a the CMS crystal
                                                       
      GPMaterial *ecal_mat =  reserveThing<GPMaterial>();
      GPMaterial_Constructor_ByElement(ecal_mat,8.28*g/cm3);

      GPElement ele_Pb;
      GPElement ele_W;
      GPElement ele_O4;

      GPElement_Constructor(&ele_Pb,82,207.2*g/mole);
      GPElement_Constructor(&ele_W,74,183.84*g/mole);
      GPElement_Constructor(&ele_O4,8*4,15.9994*4*g/mole);

      GPMaterial_AddElement(ecal_mat,ele_Pb,0.45532661);
      GPMaterial_AddElement(ecal_mat,ele_W, 0.40403397);
      GPMaterial_AddElement(ecal_mat,ele_O4,0.14063942);

      
      for ( int j = 0; j < crystal_nz ; ++j ) {
         for ( int i = 0; i < crystal_nphi ; ++i ) {
            
            iptr = i+j*crystal_nphi;
            
            GPTubs *ecal = reserveThing<GPTubs>();
            
            GPTubs_Constructor(ecal, ecal_rmin, ecal_rmax,
                               ecal_dz, ecal_sphi+i*ecal_dphi,ecal_dphi);
            
            GPThreeVector ecal_trans = GPThreeVector_create(0,0,ecal_zmin+(2.0*j+1.0)*ecal_dz);
            
            GPLogicalVolume_Constructor(ecal_log+iptr, (GPVSolid*)ecal, ecal_mat);
            GPVPhysicalVolume_Constructor(ecal_phy+iptr, idRot, ecal_trans, ecal_log+iptr);
            
            //Set mother
            GPVPhysicalVolume_SetMotherLogical(ecal_phy+iptr, world_log);
            
            addLogicalVolumePointers( ecal_log+iptr);
            addPhysicalVolumePointers( ecal_phy+iptr);
            
         }
      }
      
      //add daughter volume
      for ( int j=0; j < crystal_nz ; ++j ) {
         for ( int i=0; i < crystal_nphi ; ++i ) {
            iptr = i+j*crystal_nphi;
            addLogicalVolumeDaughter( world_log, ecal_phy+iptr);
         }
      }
      
      // Register world volume pointers for relocation
      addLogicalVolumePointers( world_log );
      addPhysicalVolumePointers( world_phy );
      
      // We *known* the VP index in order of creations, let do the same
      // things.
      geomManager->fVolumesIndex = reserveNThings<GPLogicalVolume*>(crystal_n);
      addPointer(geomManager->fVolumesIndex);
      world_log->fIndex = 0;
      geomManager->fVolumesIndex[0] = world_log;
      addPointer( geomManager->fVolumesIndex[0] );
      world_phy->fIndex = 0;
      for ( int j = 0; j < crystal_nz ; ++j ) {
         for ( int i = 0; i < crystal_nphi ; ++i ) {
            iptr = i+j*crystal_nphi;
            geomManager->fVolumesIndex[iptr+1] = ecal_log+iptr;
            addPointer( geomManager->fVolumesIndex[iptr] );
            (ecal_log+iptr)->fIndex = iptr+1;
            (ecal_phy+iptr)->fIndex = iptr; // Because we have just one physical volume per logical ... and they are all in the world logical volume.
         }
      }
   }
};

#endif
