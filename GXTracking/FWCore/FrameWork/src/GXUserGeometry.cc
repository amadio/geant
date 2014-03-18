// Basic Geometry Class based on geometry_commom.hpp created by
// Otto Seiskari, August 2010

#include "GXUserGeometry.h"
#include "GPVoxelHeader.h"

void *GXUserGeometry::reserveThingUnaligned( int thingsz )
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

void *GXUserGeometry::reserveThingAligned( int thingsz, int alignment )
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

void *GXUserGeometry::reserveThing( int thingsz )
{
   return reserveThingAligned( thingsz, default_alignment );
}

void *GXUserGeometry::addThingAligned( const void *thing,
                                      int thingsz, int alignment )
{
   void *ptr = reserveThingAligned(thingsz, alignment);
   std::memcpy( ptr, thing, thingsz );
   return ptr;
}

void *GXUserGeometry::addThingUnaligned( const void *thing, int thingsz )
{
   void *ptr = reserveThingUnaligned(thingsz);
   std::memcpy( ptr, thing, thingsz );
   return ptr;
}

void *GXUserGeometry::addThing( const void *thing, int thingsz )
{
   void *ptr = reserveThing(thingsz);
   std::memcpy( ptr, thing, thingsz );
   return ptr;
}


void GXUserGeometry::addLogicalVolumePointers( GPLogicalVolume *vol )
{
   assert( vol->fMaterial != NULL );
   assert( vol->fSolid != NULL );
   
   addPointer( vol->fMaterial );
   addPointer( vol->fSolid );
   addPointer( vol->fDaughters );
   for (int i=0; i<vol->fNoDaughters; ++i )
      addPointer( vol->fDaughters[i] );
   
   addPointer( vol->fVoxel );
}

void GXUserGeometry::addPhysicalVolumePointers( GPVPhysicalVolume *vol )
{
   assert( vol->flogical != NULL );
   
   addPointer( vol->flogical );
   addPointer( vol->flmother );
}

/**
 * Add daughter to logical volume data structure
 * the daughters for one volume must be added consecutively,
 * no other "add" or "reserve" functions may be called inbetween
 */
void GXUserGeometry::addLogicalVolumeDaughter( GPLogicalVolume *vol,
                                              GPVPhysicalVolume *d )
{
   // Constructing a C-array --> no padding / alignment!
   GPVPhysicalVolume **dptr = addThingUnaligned( d );
   
   if (vol->fDaughters == NULL)
   {
      vol->fDaughters = dptr;
   }
   vol->fNoDaughters++;
}

GPSmartVoxelNode *GXUserGeometry::storeVoxelNode( const GPSmartVoxelNode &original )
{
   nvoxelnodes++;
   
   GPSmartVoxelNode *node = addThing(original);
   if ( original.fcontents != NULL ) {
      node->fcontents = static_cast<G4int*>(addThing( original.fcontents,
                                                     sizeof(G4int)*original.fnumContents ));
      addPointer( node->fcontents );
   }
   return node;
}

GPSmartVoxelHeader *GXUserGeometry::storeVoxelHeader( const GPSmartVoxelHeader &original )
{
   GPSmartVoxelHeader *v = addThing(original);
   assert( v->fnumSlices > 0 && v->fslices != NULL );
   
   v->fslices = reserveNThingsSelfAligned<GPSmartVoxelProxy*>( v->fnumSlices );
   addPointer( v->fslices );
   
   for ( int i=0; i<v->fnumSlices; ++i ) {
      assert( original.fslices[i] != NULL );
      
      v->fslices[i] = addThing( *original.fslices[i] );
      addPointer( v->fslices[i] );
      
      if ( GPSmartVoxelProxy_IsNode( original.fslices[i] ) ) {
         v->fslices[i]->fNode = storeVoxelNode( *original.fslices[i]->fNode );
         addPointer( v->fslices[i]->fNode );
      }
      else {
         assert( original.fslices[i]->fHeader != NULL );
         
         v->fslices[i]->fHeader =
         storeVoxelHeader( *original.fslices[i]->fHeader );
         
         addPointer( v->fslices[i]->fHeader );
      }
   }
   return v;
}

void GXUserGeometry::createVoxels( GPLogicalVolume *mother, G4double smartless = kInfinity )
{
   GPSmartVoxelHeader head;
   GPSmartVoxelHeader_Constructor( &head, mother, 0, smartless );
   
   if ( head.fnumSlices > 0 ) GPLogicalVolume_SetVoxelHeader(mother, storeVoxelHeader( head ) );
   
   GPSmartVoxelHeader_Destructor( &head );
}

G4double GXUserGeometry::uniformRand()
{
   const G4int MAGN = 10000;
   return (std::rand() % MAGN) / (G4double)MAGN;
}

GPRotationMatrix GXUserGeometry::randomRot()
{
   GPRotationMatrix rot = GPRotationMatrix_create(1,0,0,
                                                  0,1,0,
                                                  0,0,1);
   const G4double angX = uniformRand() * 2.0 * M_PI;
   const G4double angY = uniformRand() * 2.0 * M_PI;
   const G4double angZ = uniformRand() * 2.0 * M_PI;
   
   GPRotationMatrix_rotateX( &rot, angX );
   GPRotationMatrix_rotateY( &rot, angY );
   GPRotationMatrix_rotateZ( &rot, angZ );
   
   return rot;
}

/**
 * Actual "create" implementation. No guards against
 * this begin called multiple times.
 */
GXUserGeometry::GXUserGeometry( void *buffer, int maxsz,
                               int def_align = sizeof(G4double) )
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

GXUserGeometry::GXUserGeometry( int maxsz )
:
userbuf(false),
create_called(false),
bufsz(0),
maxbufsz(maxsz),
default_alignment(sizeof(G4double)),
nvoxelnodes(0)
{
   if ( maxsz == 0 )
      buf = NULL;
   else
      buf = std::malloc( maxsz );
   
   assert(default_alignment > 0);
}

GXUserGeometry::~GXUserGeometry()
{
   if (!userbuf) std::free( buf );
}

void GXUserGeometry::relocate( void *newbegin )
{
   for ( unsigned i=0; i<ptrs.size(); ++i )
   {
      *(byte**)((byte*)buf+ptrs[i].destoffs) =
      (byte*)newbegin + ptrs[i].targoffs;
   }
}

void GXUserGeometry::create()
{
   if ( create_called )
      throw std::runtime_error("GXUserGeometry::create called twice");
   else create_impl();
}

template void GXUserGeometry::addPointer( GPLogicalVolume **& ptr );
