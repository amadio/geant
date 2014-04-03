//
//  GPGeomManager.h
//  
//
//  Created by Philippe Canal on 4/12/13.
//
//

#ifndef GPGeomManager_H
#define GPGeomManager_H

class GPVPhysicalVolume;
class GPLogicalVolume;

class GPGeomManager
{
public:
   GPGeomManager() : fWorldPhysical(0), fVolumesIndex(0) {}

   GPVPhysicalVolume  *fWorldPhysical;
   GPLogicalVolume   **fVolumesIndex;
};

#endif // GPGeomManager_H


