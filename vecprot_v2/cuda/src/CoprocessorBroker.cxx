#include "CoprocessorBroker.h"


// Needed by GXEMPhysicsUtils.h for now
#include <iostream>
#include <fstream>
typedef double G4double;

#include "GPConstants.h"
#include "GXFieldMap.h"
#include "GXFieldMapData.h"

#include "GXTrack.h"
#include "GPThreeVector.h"

//Geometry
#include "GPVGeometry.h"
#include "GPUserGeometry.h"
#include "GPSimpleEcal.h"
#include "GPVPhysicalVolume.h"

//EMPhysics
#include "GPConstants.h"
#include "GPPhysics2DVector.h"

#include "GPPhysicsTable.h"
#include "GPPhysicsTableType.h"

#include "GeantTrack.h"

#include "GXTrack.h"
#include "GXTrackLiason.h"

#include "GPEMPhysicsUtils.h"
//#include "GXRunManager.h"
//#include "GPEMPhysicsUtils.h"

#include <vector>
#include <math.h>

// #include <cuda_runtime.h>

#include "CoprocessorBrokerKernel.h"

#include "TGeoBranchArray.h"
#include "TGeoNode.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"

// To access the list of baskets.
#include "WorkloadManager.h"
#include "GeantScheduler.h"
#include "GeantBasket.h"

static void HandleCudaError( cudaError_t err,
                             const char *file,
                             int line ) { 
    if (err != cudaSuccess) {
       Fatal("Cuda","%s (%d) in %s at line %d\n", cudaGetErrorString( err ), err, 
             file, line );
       exit( EXIT_FAILURE );
    }   
}

#define HANDLE_CUDA_ERROR( err ) (HandleCudaError( err, __FILE__, __LINE__ ))

struct GeneralTask : public CoprocessorBroker::Task {
   GeneralTask() : Task( tracking_gpu ) {}

   const char *Name() { return "GeneralTask"; }
   bool Select(GeantTrack_v &host_track, int track)
   {
      // Currently we can only handle electron, which we pretend are the only
      // particle to have charge -1.

      return -1 == host_track.fChargeV[track];
   }
};

struct EnergyTask : public CoprocessorBroker::Task {
   double  fThresHold;
   TString fName;

   EnergyTask(double threshold) : Task( tracking_gpu ), fThresHold(threshold)
   {
      fName.Form("EnergyTask - %3g",threshold);
   }

   const char *Name() { return fName; }
   bool Select(GeantTrack_v &host_track, int track)
   {
      // Currently we can only handle electron, which we pretend are the only
      // particle to have charge -1.

      return -1 == host_track.fChargeV[track] && host_track.fEV[track] > fThresHold;
   }
};

struct EnergyElectronTask : public CoprocessorBroker::Task {
   double  fThresHold;
   TString fName;

   EnergyElectronTask(double threshold) : Task( electron_multistage_gpu ), fThresHold(threshold)
   {
      fName.Form("EnergyElectronTask - %3g",threshold);
   }

   const char *Name() { return fName; }
   bool Select(GeantTrack_v &host_track, int track)
   {
      // Currently we can only handle electron, which we pretend are the only
      // particle to have charge -1.

      return -1 == host_track.fChargeV[track] && host_track.fEV[track] > fThresHold;
   }
};

struct EnergyElectronSingleTask : public CoprocessorBroker::Task {
   double  fThresHold;
   TString fName;

   EnergyElectronSingleTask(double threshold) : Task( electron_gpu ), fThresHold(threshold)
   {
      fName.Form("EnergyElectronSingleTask - %3g",threshold);
   }

   const char *Name() { return fName; }
   bool Select(GeantTrack_v &host_track, int track)
   {
      // Currently we can only handle electron, which we pretend are the only
      // particle to have charge -1.

      return -1 == host_track.fChargeV[track] && host_track.fEV[track] > fThresHold;
   }
};


#if 0
struct ElectronTask : public CoprocessorBroker::Task {
   ElectronTask() : Task( tracking_electron_gpu ) {}

   const char *Name() { return "ElectronTask"; }
   bool Select(GeantTrack_v &host_track, int track)
   {
      // Currently we can only handle electron, which we pretend are the only
      // particle to have charge -1.

      return -1 == host_track.fChargeV[track];
   }
};

struct PhotonTask : public CoprocessorBroker::Task {
   PhotonTask() : Task( tracking_photon_gpu ) {}

   const char *Name() { return "PhotonTask"; }
   bool Select(GeantTrack_v &host_track, int track)
   {
      // Currently we can only handle photon, which we pretend are the only
      // particle to have charge 0.

      return 0 == host_track.fChargeV[track];
   }
};
#endif

void DevicePtrBase::Malloc(unsigned long size) {
   HANDLE_CUDA_ERROR( cudaMalloc((void**)&fPtr,size) );
}

DevicePtrBase::~DevicePtrBase() {
   if (fPtr) HANDLE_CUDA_ERROR( cudaFree(fPtr) );
}

void DevicePtrBase::MemcpyToDevice(const void* what, unsigned long nbytes)
{
   HANDLE_CUDA_ERROR(cudaMemcpy(fPtr,what,nbytes,
                                cudaMemcpyHostToDevice) );
}

void DevicePtrBase::MemcpyToHostAsync(void* where, unsigned long nbytes, cudaStream_t stream)
{
   HANDLE_CUDA_ERROR(cudaMemcpyAsync(where, fPtr, nbytes, cudaMemcpyDeviceToHost, stream));
}

void SecondariesTable::Alloc(size_t maxTracks)
{
   fDevStackSize.Alloc();
   //fDevOffset.Alloc();
   fDevTracks.Alloc(maxTracks * kMaxNumStep * kMaxSecondaryPerStep);
}

GXFieldMap **ReadFieldMap(const char *fieldMapFile)
{
   GXFieldMap** fieldMap;
   fieldMap = (GXFieldMap **) malloc (nbinZ*sizeof(GXFieldMap *));
   for (int j = 0 ; j < nbinZ ; j++) {
      fieldMap[j] = (GXFieldMap *) malloc (nbinR*sizeof(GXFieldMap));
   }
   std::ifstream ifile(fieldMapFile, std::ios::in | std::ios::binary | std::ios::ate);
   if (ifile.is_open()) {

      //field map structure
      GXFieldMapData fd;

      std::ifstream::pos_type fsize = ifile.tellg();
      size_t dsize = sizeof(GXFieldMapData);

      long int ngrid = fsize/dsize;
      ifile.seekg (0,  std::ios::beg);

      std::cout << "... transportation ... Loading magnetic field map: "
      << fieldMapFile << std::endl;

      for(int i = 0 ; i < ngrid ; i++) {
         ifile.read((char *)&fd, sizeof(GXFieldMapData));

         //check validity of input data
         if(abs(fd.iz) > noffZ || fd.ir > nbinR) {
            std::cout << " Field Map Array Out of Range" << std::endl;
         }
         else {
            fieldMap[noffZ+fd.iz][fd.ir].Bz = fd.Bz;
            fieldMap[noffZ+fd.iz][fd.ir].Br = fd.Br;
         }
      }
      ifile.close();
   } else {
      Error("ReadFieldMap","Could not open the file %s\n",fieldMapFile);
   }
   return fieldMap;
}

bool ReadPhysicsTable(GPPhysicsTable &table, const char *filename, bool useSpline)
{

   readTable(&table,filename);
   unsigned int nv = table.nPhysicsVector;
   for(unsigned int j=0; j < nv; j++){
      table.physicsVectors[j].SetSpline(useSpline);
   }
   //  GPPhysicsTable_Print(&table);
   return true;
}


CoprocessorBroker::TaskData::TaskData() : fTrack(0),fTrackId(0),fPhysIndex(0),fLogIndex(0),fNStaged(0),
                                          fPrioritizer(0),
                                          fStreamId(0),
                                          fThreadId(-1),
                                          fBasket(0), fQueue(0)
{
   // Default constructor.
}

CoprocessorBroker::TaskData::~TaskData() {
   // Destructor.

   // We do not own fQueue.
   // We do not own fTrackCollection.
   // We do not own the fBasket (or do we?)

   delete [] fTrack;
   delete [] fTrackId;
   delete [] fPhysIndex;
   delete [] fLogIndex;
   HANDLE_CUDA_ERROR( cudaStreamDestroy(fStream) );
}


bool CoprocessorBroker::TaskData::CudaSetup(unsigned int streamid, int nblocks, int nthreads, int maxTrackPerThread)
{

   fStreamId = streamid;
   HANDLE_CUDA_ERROR( cudaStreamCreate(&fStream) );

   HANDLE_CUDA_ERROR( cudaMemcpyToSymbol("gPropagator_fBmag", &(gPropagator->fBmag), sizeof(double), size_t(0), cudaMemcpyHostToDevice) );

   HANDLE_CUDA_ERROR( cudaMemcpyToSymbol("gTolerance", &(TGeoShape::Tolerance()), sizeof(double), size_t(0), cudaMemcpyHostToDevice) );

   //prepare random engines on the device
   fdRandStates.Alloc( nblocks*nthreads );
   curand_setup_gpu(fdRandStates, time(NULL), nblocks, nthreads);

   unsigned int maxTrackPerKernel = nblocks*nthreads*maxTrackPerThread;
   fChunkSize = maxTrackPerKernel;
   fDevTrack.Alloc(maxTrackPerKernel);
   fDevAltTrack.Alloc(maxTrackPerKernel);
   fDevTrackPhysIndex.Alloc(maxTrackPerKernel);
   fDevTrackLogIndex.Alloc(maxTrackPerKernel);
   fDevSecondaries.Alloc(maxTrackPerKernel);

   fDevScratch.Alloc(1);
   fDevTrackTempSpace.Alloc(maxTrackPerKernel);

   fTrack = new GXTrack[maxTrackPerKernel];
   fTrackId = new int[maxTrackPerKernel];
   fPhysIndex = new int[maxTrackPerKernel];
   fLogIndex = new int[maxTrackPerKernel];


   return true;
}

void CoprocessorBroker::TaskData::Push(concurrent_queue *q)
{
   // Add this helper to the queue and record it as the
   // default queue.

   if (q) {
      fQueue = q;
      fQueue->push(this);
   } else if (fQueue) {
      fQueue->push(this);
   }
}

void CoprocessorBroker::TaskData::Reset()
{
   // Get the stream helper ready for re-use.
   // by setting fNStaged to zero and add the helper
   // to the list of available helpers.

   ResetNStaged();
   Push();
}

CoprocessorBroker::CoprocessorBroker() : fdGeometry(0)
   ,fNextTaskData(0)
   ,fNblocks(0),fNthreads(0)
   ,fMaxTrackPerThread(1)
   ,fTotalWork(0)
//,fdFieldMap(0)
//,fdStates0(0),fdStates1(0),fdRandStates0(0),fdRandStates1(0)
{
   // Default constructor.

}

CoprocessorBroker::~CoprocessorBroker()
{
   HANDLE_CUDA_ERROR( cudaDeviceSynchronize() );
   for(unsigned int i = 0; i < fTasks.size(); ++i) {
      delete fTasks[i];
   }
   fTasks.clear();
   for(unsigned int i = 0; i < fTaskData.size(); ++i) {
      delete fTaskData[i];
   }
   fTaskData.clear();

   HANDLE_CUDA_ERROR( cudaFree(fdGeometry) );
//   cudaFree(fdFieldMap);
//   cudaFree(fd_eBremTable);
//   cudaFree(fd_eIoniTable);
//   cudaFree(fd_mscTable);
//   cudaFree(fdRandStates0);
//   cudaFree(fdRandStates1);

}

bool CoprocessorBroker::UploadGeometry(GPVGeometry *geom)
{
   // Prepare the geometry for the device and upload it to the device's memory.

   HANDLE_CUDA_ERROR(cudaMalloc( (void**)&fdGeometry, geom->size() ));
   geom->relocate( fdGeometry );
   HANDLE_CUDA_ERROR(cudaMemcpy(fdGeometry, geom->getBuffer(), geom->size(), cudaMemcpyHostToDevice));

   return true;
}

bool CoprocessorBroker::UploadMagneticField(GXFieldMap** fieldMap)
{

   printf("Creating magnetic field map on the GPU device\n");
   //prepare fiesldMap array: fieldMap[nbinZ][nbinR];
   GXFieldMap *bmap_h = (GXFieldMap *) malloc (nbinZ*nbinR*sizeof(GXFieldMap));
   for (int i = 0 ; i < nbinZ ; i++) {
      for (int j = 0 ; j < nbinR ; j++) {
         bmap_h[i+j*nbinZ].Bz = fieldMap[i][j].Bz;
         bmap_h[i+j*nbinZ].Br = fieldMap[i][j].Br;
      }
   }

   printf("Copying data from host to device\n");

//   cudaMalloc((void**)&fdFieldMap, nbinZ*nbinR*sizeof(GXFieldMap)) ;
//   cudaMemcpy(fdFieldMap,bmap_h,nbinZ*nbinR*sizeof(GXFieldMap),
//              cudaMemcpyHostToDevice);
   fdFieldMap.Alloc(nbinZ*nbinR);
   fdFieldMap.ToDevice(bmap_h,nbinZ*nbinR);
   free(bmap_h);

   return true;
}

bool CoprocessorBroker::UploadPhysicsTable(const GPPhysicsTable *table, unsigned int nTables, GPPhysics2DVector* sbData, size_t maxZ)
{
   fd_PhysicsTable.Alloc(nTables);
   fd_PhysicsTable.ToDevice(table,nTables);

   fd_SeltzerBergerData.Alloc(maxZ);
   fd_SeltzerBergerData.ToDevice(sbData,maxZ);

   return true;
}

void setup(CoprocessorBroker *broker,
           int nphi = 4,
           int nz   = 3,
           double density = 8.28)
{
   //1. Create the geometry.
   GPVGeometry *geom = new GPSimpleEcal(nphi,nz,density);
   geom->create();

   broker->UploadGeometry(geom);

   //2. Read magnetic field map
   const char* fieldMapFile = getenv("GP_BFIELD_MAP");
   fieldMapFile = (fieldMapFile) ? fieldMapFile : "data/cmsExp.mag.3_8T";
   GXFieldMap** fieldMap = ReadFieldMap(fieldMapFile);

   //3. Create magnetic field on the device
   broker->UploadMagneticField(fieldMap);

   // 4. Prepare EM physics tables
   GPPhysicsTable physicsTable[kNumberPhysicsTable];

   TString filename;
   for(int it = 0 ; it < kNumberPhysicsTable ; ++it) {
      filename.Form("data/%s",GPPhysicsTableName[it]);
      readTableAndSetSpline(&physicsTable[it],filename);
   }

   //G4SeltzerBergerModel data
   unsigned int maxZ = 92;
   GPPhysics2DVector* sbData = new GPPhysics2DVector[maxZ];

   char sbDataFile[256];
   for(unsigned int iZ = 0 ; iZ < maxZ ; ++iZ) {
      sprintf(sbDataFile,"data/brem_SB/br%d",iZ+1);
      std::ifstream fin(sbDataFile);
      G4bool check = RetrieveSeltzerBergerData(fin, &sbData[iZ]);
      if(!check) {
         printf("Failed To open SeltzerBerger Data file for Z= %d\n",iZ+1);
      }
   }
   printf("Copying GPPhysicsTable data from host to device\n");
   broker->UploadPhysicsTable(physicsTable,kNumberPhysicsTable,sbData,maxZ);

   delete [] sbData;

}

bool CoprocessorBroker::CudaSetup(int nblocks, int nthreads, int maxTrackPerThread)
{
   int deviceCount = 0;
   cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

   if (error_id != cudaSuccess) {
     Printf("Cuda CoprocessorBroker disabled because of issue #%d:%s\n", (int)error_id, cudaGetErrorString(error_id));
     return false;
   }

   setup(this);

   fNblocks = nblocks;
   fNthreads = nthreads;
   fMaxTrackPerThread = maxTrackPerThread;

   //fTasks.push_back(new GeneralTask());
   fTasks.push_back(new EnergyElectronTask(6));
   fTasks.push_back(new EnergyElectronTask(4));
   fTasks.push_back(new EnergyElectronTask(2));
   fTasks.push_back(new EnergyElectronTask(0));

   //initialize the stream
   for(unsigned int i=0; i < 2+fTasks.size(); ++i) {
      TaskData *data = new TaskData();
      data->CudaSetup(i,nblocks,nthreads,maxTrackPerThread);
      data->Push(&fHelpers);
      fTaskData.push_back(data);
   }
   return true;
}

// bool CanHandleTrack(GeantTrack &track)
// {
//    // Currently we can only handle electron, which we pretend are the only
//    // particle to have charge -1.

//    return -1 == track.charge;
// }

unsigned int CoprocessorBroker::TaskData::TrackToDevice(CoprocessorBroker::Task *task,
                                                        int tid,
                                                        GeantBasket &basket,
                                                        unsigned int startIdx)
{
   if (fThreadId == -1) {
      fThreadId = tid;
   } else if (fThreadId != tid) {
      Fatal("TrackToDevice","The stream %p is already assigned to thread %d when requested to do work for thread %d",
            this, fThreadId, tid);
      return 0;
   }

   if (!fBasket) {
      GeantBasketMgr *bmgr = basket.GetBasketMgr();
      fBasket = bmgr->GetNextBasket();
   }

   unsigned int start = fNStaged;
   unsigned int count = 0;
   GeantTrack_v &input = basket.GetInputTracks();
   unsigned int basketSize = input.GetNtracks();
   for(unsigned int hostIdx = startIdx;
       fNStaged < fChunkSize && hostIdx < basketSize;
       ++hostIdx ) {

      ++count;

      if (input.fHoles->TestBitNumber(hostIdx))
         continue;

      VolumePath_t *path = input.fPathV[hostIdx];
//      if (path->GetLevel()>1) {
//         fprintf(stderr,"DEBUG: for %d level is %d\n",trackin[hostIdx],path->GetLevel());
//      }
      //fprintf(stderr,"DEBUG8: See track #%d\n",trackin[hostIdx]);
      TGeoVolume *logVol = path->GetCurrentNode()->GetVolume();
      int volumeIndex = ((GeantBasketMgr*)(logVol->GetFWExtension()))->GetNumber();

      if (task->Select(input,hostIdx)) {

         fTrack[fNStaged].x  = input.fXposV[hostIdx];
         fTrack[fNStaged].y  = input.fYposV[hostIdx];
         fTrack[fNStaged].z  = input.fZposV[hostIdx];
         fTrack[fNStaged].px = input.fXdirV[hostIdx];
         fTrack[fNStaged].py = input.fYdirV[hostIdx];
         fTrack[fNStaged].pz = input.fZdirV[hostIdx];
         fTrack[fNStaged].q  = input.fChargeV[hostIdx];
         fTrack[fNStaged].s  = input.fStepV[hostIdx];
         fTrack[fNStaged].E  = input.fEV[hostIdx];
         if (fTrack[fNStaged].s == 0) {
            // humm cheat for now :(
            fTrack[fNStaged].s = 1.0+1*(2.0*rand()/RAND_MAX-1.0);
         }
         fLogIndex[fNStaged] = volumeIndex;

         // Finds the node index.
         int nodeIndex;
         if (gGeoManager->GetTopVolume() == logVol || path->GetLevel() == 0) {
            nodeIndex = 0; // humm why?
         } else {
            TGeoVolume *mother = path->GetNode(path->GetLevel()-1)->GetVolume();
            nodeIndex = mother->GetNodes()->IndexOf(path->GetCurrentNode());
         }
         fPhysIndex[fNStaged] = nodeIndex;

         fTrackId[fNStaged] = input.PostponeTrack(hostIdx,fBasket->GetOutputTracks());

         ++fNStaged;
      }
   }

   // count = min(fChunkSize,basketSize-startIdx);
   int ntrack = fNStaged - start;
   HANDLE_CUDA_ERROR(cudaMemcpyAsync(fDevTrack+start, fTrack+start, ntrack*sizeof(GXTrack),
                                     cudaMemcpyHostToDevice, fStream));
   HANDLE_CUDA_ERROR(cudaMemcpyAsync(fDevTrackLogIndex+start, fLogIndex+start, ntrack*sizeof(int),
                                     cudaMemcpyHostToDevice, fStream));
   HANDLE_CUDA_ERROR(cudaMemcpyAsync(fDevTrackPhysIndex+start, fPhysIndex+start, ntrack*sizeof(int),
                                     cudaMemcpyHostToDevice, fStream));

   return count;
}

unsigned int CoprocessorBroker::TaskData::TrackToHost()
{
   WorkloadManager *mgr = WorkloadManager::Instance();
   GeantScheduler *sch = mgr->GetScheduler();
   condition_locker &sched_locker = mgr->GetSchLocker();
   std::vector<TGeoNode *> array;
   int last_logical = -1;
   int last_phys = -1;
   GeantTrack_v &output = fBasket->GetOutputTracks();
   // unsigned int basketSize = output.GetNtracks();
   for(unsigned int devIdx = 0;
       devIdx < fNStaged;
       ++devIdx) {
      if ( fTrackId[devIdx] < 0 ) {
         // fprintf(stderr,"DEBUG: missig a track %d %d\n",devIdx,fTrackId[devIdx]);
         continue;
      }
      output.fXposV[fTrackId[devIdx]] = fTrack[devIdx].x;
      output.fYposV[fTrackId[devIdx]] = fTrack[devIdx].y;
      output.fZposV[fTrackId[devIdx]] = fTrack[devIdx].z;
      output.fXdirV[fTrackId[devIdx]] = fTrack[devIdx].px;
      output.fYdirV[fTrackId[devIdx]] = fTrack[devIdx].py;
      output.fZdirV[fTrackId[devIdx]] = fTrack[devIdx].pz;
      output.fChargeV[fTrackId[devIdx]] = fTrack[devIdx].q;
      output.fStepV[fTrackId[devIdx]] = fTrack[devIdx].s;
      output.fEV[fTrackId[devIdx]] = fTrack[devIdx].E;

      // NOTE: need to update the path ....

      if (fLogIndex[devIdx] == -1) {
         // The track/particle is done.
         // fprintf(stderr,"DEBUG: stopping %d in basket %d %d:%d - isalive %d\n", fTrackId[devIdx],fLogIndex[devIdx],
         //         fHostTracks[fTrackId[devIdx]]->event,fTrackId[devIdx],fHostTracks[fTrackId[devIdx]]->IsAlive())
         //gPropagator->StopTrack(gPropagator->fTracks[fTrackId[devIdx]]);

         // Note: could also have been kExitingSetup
         output.fStatusV[fTrackId[devIdx]] = kKilled;
      } else if (1) {
         output.fStatusV[fTrackId[devIdx]] = kAlive;

         if (array.size() && last_logical == fLogIndex[devIdx] && last_phys == fPhysIndex[devIdx]) {
            // Use the already setup array
         } else {
            // NOTE: to do a full rebuild, I would really need to calculate
            // the path (in term of indexes and matrix) on the GPU.
            last_phys = fPhysIndex[devIdx];
            last_logical = fLogIndex[devIdx];
            array.clear();

            TGeoNode *node = gGeoManager->GetTopNode();
            array.push_back( node );

            // NOTE: cheating .. what we really need to do is to find the
            // logical volume/node for each and every level.

            // For each level do something like::
            {
               TGeoVolume *volume = node->GetVolume();
               if (!volume->GetNodes()) {
                  Fatal("CoprocessorBroker::TaskData::TrackToHost","Unexpectedly in a leaf node...");
               }
               // Was:
               //TGeoVolume *sub_volume = mgr->GetBasketArray()[fLogIndex[devIdx]]->GetVolume();
               // Could be:
               //TGeoVolume *sub_volume = wm->GetSchduler()->fBasketMgr[fLogIndex[devIdx]]->GetVolume();
               TGeoVolume *sub_volume = (TGeoVolume*) gGeoManager->GetListOfVolumes()->At(fLogIndex[devIdx]);
               if (sub_volume != volume) {

                  node = (TGeoNode*) volume->GetNodes()->At( fPhysIndex[devIdx]);
                  array.push_back( node );
               }
            }
         }
#if 1 // from GPU
         // NOTE: somehow adding the navigator makes the result much more stable!
         TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
         if (!nav) nav = gGeoManager->AddNavigator();

#ifdef USE_VECGEOM_NAVIGATOR
         // Crap ...
         // output.fPathV[fTrackId[devIdx]]->Init(&(array[0]),matrix,array.size()-1);
#else
         static TGeoIdentity *matrix = new TGeoIdentity(); // NOTE: oh well ... later ...
         output.fPathV[fTrackId[devIdx]]->Init(&(array[0]),matrix,array.size()-1);
#endif
         //fHostTracks[fTrackId[devIdx]]->fPath->Init(&(array[0]),&matrix,array.size()-1);
#else // from CPU
         TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
         if (!nav) nav = gGeoManager->AddNavigator();
         TGeoNode *f_node = nav->FindNode(fTrack[devIdx].x, fTrack[devIdx].y, fTrack[devIdx].z);
         fHostTracks[fTrackId[devIdx]]->fPath->InitFromNavigator(gGeoManager->GetCurrentNavigator());
#endif
         VolumePath_t *path = output.fPathV[fTrackId[devIdx]];
            //fHostTracks[fTrackId[devIdx]]->fPath;
//         if (fHostTracks[fTrackId[devIdx]]->fFrombdr)
//            path = fHostTracks[fTrackId[devIdx]]->fNextPath;

         if (path->GetLevel() != (long)(array.size())-1) {
            printf("Problem (bdr=%d) with the level %d vs %d\n",
                   output.fFrombdrV[fTrackId[devIdx]],path->GetLevel(),(int)array.size()-1);
         }
#ifdef USE_VECGEOM_NAVIGATOR
         const TGeoNode *pnode = path->GetNode(0);
         if (pnode != array[0]) {
            printf("Problem (bdr=%d) with the first entry %p vs %p level = %d\n",
                   output.fFrombdrV[fTrackId[devIdx]], pnode, array[0], path->GetLevel());
            pnode->Print();
            array[0]->Print();
         }
         
         if (path->GetLevel() >= 1 && path->GetNode(1) != array[1]) {
            pnode = path->GetNode(0);
            printf("Problem (bdr=%d) with the second entry %p vs %p level = %d\n", output.fFrombdrV[fTrackId[devIdx]],pnode, array[1],path->GetLevel());
            pnode->Print();
            array[1]->Print();
         }
#else
         if (path->GetArray()[0] != array[0]) {
            printf("Problem (bdr=%d) with the first entry %p vs %p level = %d\n",
                   output.fFrombdrV[fTrackId[devIdx]], path->GetArray()[0], array[0], path->GetLevel());
            path->GetArray()[0]->Print();
            array[0]->Print();
         }
         if (path->GetLevel() >= 1 && path->GetArray()[1] != array[1]) {
            printf("Problem (bdr=%d) with the second entry %p vs %p level = %d\n", output.fFrombdrV[fTrackId[devIdx]],path->GetArray()[1], array[1],path->GetLevel());
            path->GetArray()[1]->Print();
            array[1]->Print();
         }
#endif


         // Let's reschedule it.
         // fprintf(stderr,"DEBUG: rescheduling %d in basket %d\n",fTrackId[devIdx],fLogIndex[devIdx]);
         // fTrackCollection->AddTrack(fTrackId[devIdx],
                                    // mgr->GetBasketArray()[fLogIndex[devIdx]]);
      }
   }

   Int_t ntot = 0;
   Int_t nnew = 0;
   Int_t nkilled = 0;
   /* Int_t ninjected = */ sch->AddTracks(fBasket, ntot, nnew, nkilled, fPrioritizer);
   (void)ntot;
   (void)nnew;
   (void)nkilled;
   //mgr->TransportedQueue()->push(fBasket);
   sched_locker.StartOne();
   fBasket->Recycle();
   fBasket = 0;
   fThreadId = -1;
   return fNStaged;
}

// typedef void(CUDART_CB * cudaStreamCallback_t)(cudaStream_t stream, cudaError_t status, void
//                                                *userData)

void StreamReset(cudaStream_t /* stream */, cudaError_t status, void *userData)
{
   CoprocessorBroker::TaskData *helper = (CoprocessorBroker::TaskData*)userData;
   helper->Reset();
}

void TrackToHost(cudaStream_t /* stream */, cudaError_t status, void *userData)
{
   CoprocessorBroker::TaskData *helper = (CoprocessorBroker::TaskData*)userData;
   helper->TrackToHost();
}


CoprocessorBroker::Stream CoprocessorBroker::GetNextStream()
{
   // Return the current stream (one that we can still add input to)
   // or return new one.

   if (!fNextTaskData) {
      fNextTaskData = dynamic_cast<TaskData*>(fHelpers.wait_and_pop());
      if (!fNextTaskData) {
         // nothing we can do at the moment
         return 0;
      }
   }
   return fNextTaskData;
}

CoprocessorBroker::Stream CoprocessorBroker::launchTask(Task *task, bool wait /* = false */)
{

   Stream stream = task->fCurrent;
   task->fIdles = 0;
   task->fCycles = 0;
   task->fCurrent = 0;

   //Printf("(%d - GPU) == Starting kernel for task %s using stream %d with %d tracks\n",
   //       stream->fThreadId, task->Name(), stream->fStreamId, stream->fNStaged );

   fTotalWork += stream->fNStaged;
   int result = task->fKernel(stream->fdRandStates,
                              1 /* nSteps */,
                              stream->fNStaged,
                              stream->fDevTrack,
                              stream->fDevAltTrack,
                              stream->fDevTrackLogIndex,
                              stream->fDevTrackPhysIndex,
                              stream->fDevSecondaries.fDevTracks,
                              stream->fDevSecondaries.fDevStackSize,

                              &( (*stream->fDevScratch)[0] ),
                              stream->fDevTrackTempSpace,

                              (GPGeomManager*)fdGeometry, fdFieldMap,
                              fd_PhysicsTable, fd_SeltzerBergerData,

                              fNblocks,fNthreads,*stream);

   // Bring back the number of secondaries created.
   int stackSize;
   stream->fDevSecondaries.fDevStackSize.FromDevice(&stackSize, *stream);
   // stream->fDevSecondaries.fTrack.FromDevice( stream->fSecondaries, stackSize, *stream);

   // Bring back the modified tracks.
   if (result == 0) {
      stream->fDevTrack.FromDevice( stream->fTrack, stream->fNStaged, *stream );
   } else {
      stream->fDevAltTrack.FromDevice( stream->fTrack, stream->fNStaged, *stream );
   }
   stream->fDevTrackPhysIndex.FromDevice( stream->fPhysIndex, stream->fNStaged, *stream );
   stream->fDevTrackLogIndex.FromDevice( stream->fLogIndex, stream->fNStaged, *stream );

   HANDLE_CUDA_ERROR(cudaStreamAddCallback(stream->fStream, TrackToHost, stream, 0 ));
   HANDLE_CUDA_ERROR(cudaStreamAddCallback(stream->fStream, StreamReset, stream, 0 ));

   if (wait) {
      // Use this when you need to insure the printf are actually printed.
      HANDLE_CUDA_ERROR(cudaStreamSynchronize(*stream));
   }
   // cudaDeviceSynchronize

   return stream;
}

CoprocessorBroker::Stream CoprocessorBroker::launchTask(bool wait /* = false */)
{
   Stream stream = 0;
   TaskColl_t::iterator task = fTasks.begin();
   TaskColl_t::iterator end = fTasks.end();
   while( task != end ) {

      if ((*task)->fCurrent != 0 && (*task)->fCurrent->fNStaged != 0) {

         stream = launchTask((*task),wait);
      }
      ++task;
   }
   // Return the last stream used ....
   return stream;
}

void CoprocessorBroker::runTask(int threadid, GeantBasket &basket)
                                // unsigned int nTracks, int volumeIndex, GeantTrack **tracks, int *trackin)
{
   bool force = false;

   unsigned int nTracks = basket.GetNinput();
   //<<<---------------------------------------------------------->>>
//         fprintf(stderr,"round = %d count1 = %d fNchunk=%d nTracks=%d\n",
//                 round,count1,fNchunk,nTracks);
//         cudaMemcpyAsync(track_d0, track_h+i, count1*sizeof(GXTrack),
//                         cudaMemcpyHostToDevice, stream0);
//

   unsigned int trackUsed = 0;
   TaskColl_t::iterator task = fTasks.begin();
   TaskColl_t::iterator end = fTasks.end();
   while( task != end ) {
      unsigned int trackLeft  = nTracks;
      unsigned int trackStart = 0;
      while (trackLeft) {
         if ((*task)->fCurrent == 0) {
            if (fNextTaskData) {
               (*task)->fCurrent = fNextTaskData;
            } else {
               // If we do not yet have a TaskData, wait for one.
               // Consequence: if we are very busy we hang on to this
               // basket until one of the task finishes.
               (*task)->fCurrent = GetNextStream();
            }
            fNextTaskData = 0;
            if (!(*task)->fCurrent) break;
         }
         TaskData *stream = (*task)->fCurrent;

         unsigned int before = stream->fNStaged;
         unsigned int count = stream->TrackToDevice(*task,
                                                    threadid,
                                                    basket,
                                                    trackStart);
         trackUsed += stream->fNStaged-before;
         unsigned int rejected = nTracks-trackStart - (stream->fNStaged-before);

         //if (((*task)->fCycles % 10000) == 1) 
         //   Printf("(%d - GPU) ================= Task %s Stream %d Tracks: %d seen %d skipped %d accumulated %d idles %d cycles", threadid, (*task)->Name(), stream->fStreamId, count, rejected, stream->fNStaged, (*task)->fIdles, (*task)->fCycles);

         if (stream->fNStaged < stream->fChunkSize
             && !force) {
            // We do not have enough tracks

            if ( before == stream->fNStaged ) {
               ++( (*task)->fIdles );
            } else {
               (*task)->fIdles = 0;
            }
            ++( (*task)->fCycles );

            unsigned int idle = (*task)->fIdles;
            unsigned int cycle = (*task)->fCycles;
            if (stream->fNStaged                // There is something
                && idle > (fTasks.size()-1)     // We are beyond the expected/normal number of idles cycles
                && cycle > (stream->fChunkSize / nTracks)  // To fill we need at least that many cycles
                && 2*idle > cycle               // Our input rate has drop in half
                )
            {
               // Printf("(%d - GPU) ================= Launching idle Task %s Stream %d Idle=%d cycle=%d accumulated=%d", threadid, (*task)->Name(), stream->fStreamId, idle, cycle, stream->fNStaged);
               // if we did not make any progress in a while, assume there is no 'interesting' track left and schedule the kernel.
            } else {
               // Continue to wait for more data ...
               break;
            }
         }

         if (stream->fNStaged) {
            launchTask(*task);
         }
         trackLeft  -= count;
         trackStart += count;
      }
      ++task;
   }

   if (trackUsed == 0) {
      // if we did not make any progress in a while, assume there is no 'interesting' track left
      // and schedule the most loaded task.

      Task *heavy = 0;
      TaskColl_t::iterator task = fTasks.begin();
      TaskColl_t::iterator end = fTasks.end();
      while( task != end ) {
         if ((*task)->fCurrent && (*task)->fCurrent->fNStaged) {
            if (heavy == 0 || (*task)->fCurrent->fNStaged > heavy->fCurrent->fNStaged) {
               heavy = *task;
            }
         }
         ++task;
      }
      if (heavy) {
         // Printf("(%d - GPU) ================= Launching heavy Task %s Stream %d Idle=%d cycle=%d accumulated=%d", threadid, heavy->Name(), heavy->fCurrent->fStreamId, heavy->fIdles, heavy->fCycles, heavy->fCurrent->fNStaged);
         launchTask(heavy);
      }
   }

   if (trackUsed != nTracks) {
      GeantTrack_v &input  = basket.GetInputTracks();
      GeantTrack_v &output = basket.GetOutputTracks();
      for(unsigned int hostIdx = 0 ;hostIdx < nTracks;
          ++hostIdx ) {

         if (input.fHoles->TestBitNumber(hostIdx))
            continue;

         input.PostponeTrack(hostIdx,output);
      }
   }
   if (!basket.GetInputTracks().IsCompact()) basket.GetInputTracks().Compact();
}

void CoprocessorBroker::waitForTasks()
{
   // Make sure all the tasks are finished.

   for(unsigned int i=0; i < fTaskData.size(); ++i) {
      if (fTaskData[i]->fNStaged) {
         HANDLE_CUDA_ERROR(cudaStreamSynchronize(*fTaskData[i]));
      }
   }
}
