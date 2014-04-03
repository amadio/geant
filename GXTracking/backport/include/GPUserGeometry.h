// Basic Geometry Class based on geometry_commom.hpp created by
// Otto Seiskari, August 2010

#ifndef GPUSERGEOMETRY_H
#define GPUSERGEOMETRY_H

#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <assert.h>

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

#include "GPVGeometry.h"

class GPUserGeometry : public GPVGeometry
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
	
  void *reserveThingUnaligned( int thingsz )
  {
    assert(thingsz > 0);
		
    const int needed = bufsz + thingsz;
    if (needed > maxbufsz)
      {
	if (userbuf)
	  throw std::runtime_error( "Buffer size exceeded" );
				
	maxbufsz = needed;
	void *newbuf = std::realloc( buf, maxbufsz );
	if ( buf != NULL && newbuf != buf )
	  {
	    buf = newbuf;
	    throw std::runtime_error( "Buffer reallocated to different location" );
	  }
			
	buf = newbuf;
      }
    void * const bufpos = (void*)((byte*)buf+bufsz);
    bufsz = needed;
    return bufpos;
  }
	
  void *reserveThingAligned( int thingsz, int alignment )
  {
    assert(thingsz>0 && alignment>0);
		
    int padding = 0;
    if ( bufsz % alignment > 0 )
      {
	padding = alignment - (bufsz%alignment);
      }
		
    byte *m = (byte*)reserveThingUnaligned( thingsz + padding );
		
    const byte PADDING_BYTE = 0x0;
		
    // padding contents
    while ( padding-- > 0 ) *m++ = PADDING_BYTE;
    return (void*)m;
  }
	
  void *reserveThing( int thingsz )
  {
    return reserveThingAligned( thingsz, default_alignment );
  }
	
  void *addThingAligned( const void *thing, int thingsz, int alignment )
  {
    void *ptr = reserveThingAligned(thingsz, alignment);
    std::memcpy( ptr, thing, thingsz );
    return ptr;
  }
	
  void *addThingUnaligned( const void *thing, int thingsz )
  {
    void *ptr = reserveThingUnaligned(thingsz);
    std::memcpy( ptr, thing, thingsz );
    return ptr;
  }
	
  void *addThing( const void *thing, int thingsz )
  {
    void *ptr = reserveThing(thingsz);
    std::memcpy( ptr, thing, thingsz );
    return ptr;
  }
	
 protected:
	
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
	
  void addLogicalVolumePointers( G4LogicalVolume * /*vol*/ )
  {
     /* assert( vol->GetMaterial() != NULL ); */
     /* assert( vol->GetSolid() != NULL ); */
		
     /* addPointer( vol->GetMaterial() ); */
     /* addPointer( vol->GetSolid() ); */
     /* addPointer( vol->GetDaughters() ); */
     /* for (int i=0; i<vol->GetNoDaughters(); ++i ) */
     /* addPointer( vol->fDaughters[i] ); */
  }
	
  void addPhysicalVolumePointers( G4VPhysicalVolume * /*vol*/ )
  {
    /* assert( vol->flogical != NULL ); */
		
    /* addPointer( vol->flogical ); */
    /* addPointer( vol->flmother ); */
  }
	
  /**
   * Add daughter to logical volume data structure
   * the daughters for one volume must be added consecutively,
   * no other "add" or "reserve" functions may be called inbetween
   */
  void addLogicalVolumeDaughter( G4LogicalVolume * /* vol */, G4VPhysicalVolume * /*d*/ )
  {
    // Constructing a C-array --> no padding / alignment!
    /* G4VPhysicalVolume **dptr = addThingUnaligned( d ); */
		
    /* if (vol->fDaughters == NULL) */
    /*   { */
    /*     vol->fDaughters = dptr; */
    /*   } */
    /* vol->fNoDaughters++; */
  }
	
  G4double uniformRand()
  {
    const G4int MAGN = 10000;
    return (std::rand() % MAGN) / (G4double)MAGN;
  }
	
  G4RotationMatrix randomRot()
  {
     G4RotationMatrix rot; /* = G4RotationMatrix_create(1,0,0, */
			   /*      		   0,1,0, */
			   /*      		   0,0,1); */
    const G4double angX = uniformRand() * 2.0 * M_PI;
    const G4double angY = uniformRand() * 2.0 * M_PI;
    const G4double angZ = uniformRand() * 2.0 * M_PI;
		
    rot.rotateX( angX );
    rot.rotateY( angY );
    rot.rotateZ( angZ );
		
    return rot;
  }
	
  /**
   * Actual "create" implementation. No guards against
   * this begin called multiple times.
   */
  virtual void create_impl() = 0;
	
public:

  /**
   * Initialize geometry object, use the given buffer for the
   * geometry. Ownership of the buffer is not transferred.
   * 
   * @param buffer buffer to use for the geometry
   * @param maxsz maximum size of the buffer
   * @param def_align default memory alignment for geometry objects
   */
  GPUserGeometry( void *buffer, int maxsz, int def_align = sizeof(G4double) )
   :
  userbuf(true),
    create_called(false),
    buf(buffer),
    bufsz(0),
    maxbufsz(maxsz),
    default_alignment(def_align),
    nvoxelnodes(0)
      {
	assert(default_alignment > 0);
      }
	
  /**
   * Initialize geometry object and allocate a geometry buffer of
   * given size.
   * 
   * @param maxsz maximum size of the buffer
   * @param def_align default memory alignment for geometry objects
   */
  GPUserGeometry( int maxsz = 0, int def_align = sizeof(G4double) )
   :
  userbuf(false),
    create_called(false),
    bufsz(0),
    maxbufsz(maxsz),
    default_alignment(def_align),
    nvoxelnodes(0)
      {
	if ( maxsz == 0 )
	  buf = NULL;
	else
	  buf = std::malloc( maxsz );
			
	assert(default_alignment > 0);
      }
	
  ~GPUserGeometry()
    {
      if (!userbuf) std::free( buf );
    }
	
  int size() const { return bufsz; }
  void *getBuffer() { return buf; }
	
  /** Get default memory alignment for the geometry objects */
  int getAlignment() const { return default_alignment; }
	
  /** @return the total number of voxel nodes in the buffer */
  int getNumVoxelNodes() const { return nvoxelnodes; }
	
  void relocate( void *newbegin )
  {
    for ( unsigned i=0; i<ptrs.size(); ++i )
      {
	*(byte**)((byte*)buf+ptrs[i].destoffs) =
	  (byte*)newbegin + ptrs[i].targoffs;
      }
  }
	
  void create()
  {
    if ( create_called )
      throw std::runtime_error("GPUserGeometry::create called twice");
    else create_impl();
  }
};

#endif
