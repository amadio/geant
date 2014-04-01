// Geometry Base Class based on geometry.hpp created by
// Otto Seiskari, August 2010

#ifndef GPVGEOMETRY_H
#define GPVGEOMETRY_H

class GPVGeometry
{
public:
	
  typedef unsigned char byte;
	
  virtual ~GPVGeometry() {}

  // Create the geometry. Must be called exactly once.
  virtual void create() = 0;
	
  // Change the geometry buffer base address. Does not move the array but 
  // changes the internal pointers in the geometry data structure so that 
  // the buffer may be copied to a memory segment beginning from the given 
  // address.
  virtual void relocate( void *newbegin ) = 0;
	
  // Get geometry buffer size in bytes
  virtual int size() const = 0;
	
  // Get a pointer to the geometry buffer
  virtual void *getBuffer() = 0;
	
  // Get the scale (some dimension) of the geometry
  virtual double getScale() const = 0;

  // Get the total number of voxel nodes in the buffer
  virtual int getNumVoxelNodes() const { return 0; }
	
};

#endif
