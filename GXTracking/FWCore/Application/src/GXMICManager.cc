#include "GXMICManager.h"
#include "malloc.h"

GXMICManager::GXMICManager() : GXVCoprocessorManager("MIC")
{
  ;
}

GXMICManager::~GXMICManager()
{
  ;
}

void GXMICManager::Initialize()
{
  // if offload compilation is enabled
  int num_devices = 0;
  int device_num = 0;

#ifdef __INTEL_OFFLOAD
  printf("--- Intel(R) Xeon Phi(TM) Devices ---\n");
  num_devices = _Offload_number_of_devices();
  printf("--- Number of Target Devices: %d\n",num_devices);

  device_num = _Offload_get_device_number();
  printf("--- Which Device number : %d\n",device_num);
#endif

}

void GXMICManager::AllocateDeviceMemory(/* GXTaskData* taskData_h */
					GXVGeometry *geom,
					GXFieldMap** fieldMap2D,
					GXPhysicsTable* physicsTable,
					GXPhysics2DVector* sbData)
{
  //Magnetic Field Map
  fieldMap_d = (GXFieldMap *) _mm_malloc(nbinZ*nbinR*sizeof(GXFieldMap),128);

  for (int i = 0 ; i < nbinZ ; i++) {
    for (int j = 0 ; j < nbinR ; j++) {
      fieldMap_d[i+j*nbinZ].Bz = fieldMap2D[i][j].Bz;
      fieldMap_d[i+j*nbinZ].Br = fieldMap2D[i][j].Br;
    }
  }

  geom_d = (GPVGeometry::byte*) _mm_malloc (geom->size(),128) ;
  geom->relocate( geom_d );
  memcpy(geom_d,geom->getBuffer(),geom->size());  

 //physics tables
  physicsTable_d =
   (GXPhysicsTable*)_mm_malloc(kNumberPhysicsTable*sizeof(GXPhysicsTable),128);
  memcpy(physicsTable_d,physicsTable,
	 kNumberPhysicsTable*sizeof(GXPhysicsTable));
 
  //DB data
  sbData_d = 
    (GXPhysics2DVector*) _mm_malloc(maxElements*sizeof(GXPhysics2DVector),128);
  memcpy(sbData_d,sbData,sizeof(maxElements*sizeof(GXPhysics2DVector)));

}

void GXMICManager::DeallocateDeviceMemory()
{

}

void GXMICManager::UploadTaskData()
{
  ;
}
