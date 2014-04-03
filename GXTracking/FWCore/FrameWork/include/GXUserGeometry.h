// Basic Geometry Class based on geometry_commom.hpp created by
// Otto Seiskari, August 2010

#ifndef GXUSERGEOMETRY_H
#define GXUSERGEOMETRY_H

#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <assert.h>

#include "GPLogicalVolume.h"
#include "GPVPhysicalVolume.h"
#include "GPSmartVoxelHeader.h"
#include "GPSmartVoxelNode.h"
#include "GPSmartVoxelProxy.h"

#include "GXVGeometry.h"

class GXUserGeometry : public GXVGeometry
{
 private:
  typedef struct
  {
    int destoffs;
    int targoffs;
  }
  StrayPointer;

  const bool userbuf;
  bool create_called;
	
  void *buf;
  int bufsz, maxbufsz;
  int default_alignment;
	
  int nvoxelnodes;
	
  std::vector< StrayPointer > ptrs;
	
  void *reserveThingUnaligned( int thingsz );
	
  void *reserveThingAligned( int thingsz, int alignment );
	
  void *reserveThing( int thingsz );
	
  void *addThingAligned( const void *thing, int thingsz, int alignment );
  void *addThingUnaligned( const void *thing, int thingsz );
	
  void *addThing( const void *thing, int thingsz );
	
 protected:
	
  template <class T> void addPointer( T *& ptr )
  {
     if ( ptr != NULL )
	{
           byte **p = (byte**)(&ptr);
           StrayPointer s = {
              (byte*)p - (byte*)buf,
              *p - (byte*)buf
           };
           ptrs.push_back(s);
	}
  }

  template <class T> T* reserveThing()
    {
      return static_cast<T*>(reserveThing(sizeof(T)));
    }
	
  template <class T> T* reserveThingAligned( int alignment = sizeof(T) )
    {
      return static_cast<T*>(reserveThingAligned(sizeof(T), alignment));
    }
	
  template <class T> T* reserveThingUnaligned()
    {
      return static_cast<T*>(reserveThingUnaligned(sizeof(T)));
    }
	
  template <class T> T* reserveNThingsUnaligned( int N )
    {
      return static_cast<T*>(reserveThingUnaligned(sizeof(T)*N));
    }
	
  template <class T> T* reserveNThingsAligned( int N, int alignment = sizeof(T) )
    {
      if ( sizeof(T) % alignment )
	throw std::runtime_error( "Incorrect alignment for array" );
			
      return static_cast<T*>(reserveThingAligned(sizeof(T)*N, alignment));
    }
	
  template <class T> T* reserveNThingsSelfAligned( int N )
    {
      return reserveNThingsAligned<T>( N );
    }
	
  template <class T> T* reserveNThings( int N )
    {
      return reserveNThingsAligned<T>( N, default_alignment );
    }
	
  template <class T> T* addThing( const T &thing )
    {
      return static_cast<T*>( addThing( &thing, sizeof(T) ) );
    }
	
  template <class T> T* addThingAligned( const T &thing, int alignment = sizeof(T) )
    {
      return static_cast<T*>( addThingAligned( &thing, sizeof(T), alignment ) );
    }
	
  template <class T> T* addThingUnaligned( const T &thing )
    {
      return static_cast<T*>( addThingUnaligned( &thing, sizeof(T) ) );
    }
	
  void addLogicalVolumePointers( GPLogicalVolume *vol );
	
  void addPhysicalVolumePointers( GPVPhysicalVolume *vol );
	
  void addLogicalVolumeDaughter( GPLogicalVolume *vol, GPVPhysicalVolume *d );

  GPSmartVoxelNode* storeVoxelNode( const GPSmartVoxelNode &original );
  GPSmartVoxelHeader* storeVoxelHeader( const GPSmartVoxelHeader &original );
  void createVoxels( GPLogicalVolume *mother, G4double smartless);

  G4double uniformRand();
  GPRotationMatrix randomRot();
	
  /**
   * Actual "create" implementation. No guards against
   * this begin called multiple times.
   */
  virtual void create_impl() = 0;
	
public:

  GXUserGeometry( void *buffer, int maxsz, int def_align);
  //  GXUserGeometry( void *buffer, int maxsz, int def_align = sizeof(G4double) );
  //  GXUserGeometry( int maxsz, int def_align );
  GXUserGeometry( int maxsz );
  ~GXUserGeometry();

  int size() const { return bufsz; }
  void *getBuffer() { return buf; }
	
  /** Get default memory alignment for the geometry objects */
  int getAlignment() const { return default_alignment; }
	
  /** @return the total number of voxel nodes in the buffer */
  int getNumVoxelNodes() const { return nvoxelnodes; }
	
  void relocate( void *newbegin );
  void create();
};

#endif
