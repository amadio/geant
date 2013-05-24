#include "CoprocessorBroker.h"


// Needed by GXEMPhysicsUtils.h for now
#include <iostream>
#include <fstream>
typedef double G4double;

#include "GXEMPhysicsUtils.h"


#include "GPConstants.h"
#include "GPFieldMap.h"
#include "GPFieldMapData.h"

#include "GXTrack.h"
#include "GPThreeVector.h"

//Geometry
#include "GPVGeometry.h"
#include "GPUserGeometry.h"
#include "GPSimpleEcal.h"
#include "GPVPhysicalVolume.h"

//EMPhysics
#include "GPConstants.h"
#include "GXEMPhysicsUtils.h"

#include "GPPhysicsTable.h"

#include "GeantTrack.h"

#include "GXTrack.h"

#include <vector>
#include <math.h>

// #include <cuda_runtime.h>

#include "CoprocessorBrokerKernel.h"

#include "GeantVolumeBasket.h"
#include "TGeoBranchArray.h"
#include "TGeoNode.h"
#include "TGeoManager.h"

// To access the list of baskets.
#include "WorkloadManager.h"

void DevicePtrBase::Malloc(unsigned long size) {
   cudaMalloc((void**)&fPtr,size);
}

DevicePtrBase::~DevicePtrBase() {
   if (fPtr) cudaFree(fPtr);
}

void DevicePtrBase::MemcpyToDevice(const void* what, unsigned long nbytes)
{
   cudaMemcpy(fPtr,what,nbytes,
              cudaMemcpyHostToDevice);
}

GPFieldMap **ReadFieldMap(const char *fieldMapFile)
{
   GPFieldMap** fieldMap;
   fieldMap = (GPFieldMap **) malloc (nbinZ*sizeof(GPFieldMap *));
   for (int j = 0 ; j < nbinZ ; j++) {
      fieldMap[j] = (GPFieldMap *) malloc (nbinR*sizeof(GPFieldMap));
   }
   std::ifstream ifile(fieldMapFile, std::ios::in | std::ios::binary | std::ios::ate);
   if (ifile.is_open()) {
      
      //field map structure
      GPFieldMapData fd;
      
      std::ifstream::pos_type fsize = ifile.tellg();
      size_t dsize = sizeof(GPFieldMapData);
      
      long int ngrid = fsize/dsize;
      ifile.seekg (0,  std::ios::beg);
      
      std::cout << "... transportation ... Loading magnetic field map: "
      << fieldMapFile << std::endl;
      
      for(int i = 0 ; i < ngrid ; i++) {
         ifile.read((char *)&fd, sizeof(GPFieldMapData));
         
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


CoprocessorBroker::StreamHelper::StreamHelper() : fTrack(0),fTrackId(0),fPhysIndex(0),fLogIndex(0),fNStaged(0),
                                                  fStreamId(0),
                                                  fThreadId(-1),fHostTracks(0),
                                                  fTrackCollection(0), fQueue(0)
{
   // Default constructor.
}

CoprocessorBroker::StreamHelper::~StreamHelper() {
   // Destructor.
   
   // We do not own fHosttracks.
   // We do not own fQueue.
   // We do not own fTrackCollection.

   delete [] fTrack;
   delete [] fTrackId;
   delete [] fPhysIndex;
   delete [] fLogIndex;
   cudaStreamDestroy(fStream);
}


bool CoprocessorBroker::StreamHelper::CudaSetup(unsigned int streamid, int nblocks, int nthreads, int maxTrackPerThread)
{
   fStreamId = streamid;
   cudaStreamCreate(&fStream);
   
   //prepare random engines on the device
   fdRandStates.Alloc( nblocks*nthreads );
   curand_setup_gpu(fdRandStates, time(NULL), nblocks, nthreads);
   
   unsigned int maxTrackPerKernel = nblocks*nthreads*maxTrackPerThread;
   fChunkSize = maxTrackPerKernel;
   fDevTrack.Alloc(maxTrackPerKernel);
   fDevTrackPhysIndex.Alloc(maxTrackPerKernel);
   fDevTrackLogIndex.Alloc(maxTrackPerKernel);
   fTrack = new GXTrack[maxTrackPerKernel];
   fTrackId = new int[maxTrackPerKernel];
   fPhysIndex = new int[maxTrackPerKernel];
   fLogIndex = new int[maxTrackPerKernel];
   
   if (fTrackCollection == 0) {
      fTrackCollection = GetNextTrackCollection();
   }
   return true;
}

//______________________________________________________________________________
GeantTrackCollection *CoprocessorBroker::StreamHelper::GetNextTrackCollection()
{
   // Get a new GeantTrackCollection from the workload manager.

   GeantTrackCollection *newcoll;
   WorkloadManager *mgr = WorkloadManager::Instance();
   if (mgr->CollectorEmptyQueue()->empty()) newcoll = new GeantTrackCollection(100);
   else newcoll = (GeantTrackCollection*)mgr->CollectorEmptyQueue()->wait_and_pop();
   return newcoll;
}

void CoprocessorBroker::StreamHelper::Push(concurrent_queue *q)
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

void CoprocessorBroker::StreamHelper::Reset()
{
   // Get the stream helper ready for re-use.
   // by setting fNStaged to zero and add the helper
   // to the list of available helpers.

   ResetNStaged();
   Push();
}

CoprocessorBroker::CoprocessorBroker() : fdGeometry(0)
   ,fCurrentHelper(0)
   ,fNblocks(0),fNthreads(0)
   ,fMaxTrackPerThread(1)
   ,fKernelType(1)
   ,fTotalWork(0)
//,fdFieldMap(0)
//,fdStates0(0),fdStates1(0),fdRandStates0(0),fdRandStates1(0)
{
   // Default constructor.
   
}

CoprocessorBroker::~CoprocessorBroker()
{

   cudaFree(fdGeometry);
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
   
   cudaMalloc( (void**)&fdGeometry, geom->size() );
   geom->relocate( fdGeometry );
   cudaMemcpy(fdGeometry, geom->getBuffer(), geom->size(), cudaMemcpyHostToDevice);
   
   return true;
}

bool CoprocessorBroker::UploadMagneticField(GPFieldMap** fieldMap)
{

   printf("Creating magnetic field map on the GPU device\n");
   //prepare fiesldMap array: fieldMap[nbinZ][nbinR];
   GPFieldMap *bmap_h = (GPFieldMap *) malloc (nbinZ*nbinR*sizeof(GPFieldMap));
   for (int i = 0 ; i < nbinZ ; i++) {
      for (int j = 0 ; j < nbinR ; j++) {
         bmap_h[i+j*nbinZ].Bz = fieldMap[i][j].Bz;
         bmap_h[i+j*nbinZ].Br = fieldMap[i][j].Br;
      }
   }
   
   printf("Copying data from host to device\n");
   
//   cudaMalloc((void**)&fdFieldMap, nbinZ*nbinR*sizeof(GPFieldMap)) ;
//   cudaMemcpy(fdFieldMap,bmap_h,nbinZ*nbinR*sizeof(GPFieldMap),
//              cudaMemcpyHostToDevice);
   fdFieldMap.Alloc(nbinZ*nbinR);
   fdFieldMap.ToDevice(bmap_h,nbinZ*nbinR);
   free(bmap_h);
   
   return true;
}

bool CoprocessorBroker::Upload_eBremTable(const GPPhysicsTable &eBrem_table)
{
//   cudaMalloc((void**)&fd_eBremTable, sizeof(GPPhysicsTable));
//   cudaMemcpy(fd_eBremTable, &eBrem_table, sizeof(GPPhysicsTable),cudaMemcpyHostToDevice);
   fd_eBremTable.Alloc();
   fd_eBremTable.ToDevice(&eBrem_table);
   return true;
}

bool CoprocessorBroker::Upload_eIoniTable(const GPPhysicsTable &eIoni_table)
{
//   cudaMalloc((void**)&fd_eIoniTable, sizeof(GPPhysicsTable));
//   cudaMemcpy(fd_eIoniTable, &eIoni_table, sizeof(GPPhysicsTable),cudaMemcpyHostToDevice);
   fd_eIoniTable.Alloc();
   fd_eIoniTable.ToDevice(&eIoni_table);
   return true;
}

bool CoprocessorBroker::UploadMscTable(const GPPhysicsTable &msc_table)
{
//   cudaMalloc((void**)&fd_mscTable, sizeof(GPPhysicsTable));
//   cudaMemcpy(fd_mscTable, &msc_table, sizeof(GPPhysicsTable),cudaMemcpyHostToDevice);
   fd_mscTable.Alloc();
   fd_mscTable.ToDevice(&msc_table);
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
   fieldMapFile = (fieldMapFile) ? fieldMapFile : "cmsExp.mag.3_8T";
   GPFieldMap** fieldMap = ReadFieldMap(fieldMapFile);
   
   //3. Create magnetic field on the device
   broker->UploadMagneticField(fieldMap);

   // 4. Prepare EM physics tables
   GPPhysicsTable eBrem_table;
   ReadPhysicsTable(eBrem_table,"data/Lambda.eBrem.e-.asc",true);
   
   GPPhysicsTable eIoni_table;
   ReadPhysicsTable(eIoni_table,"data/Lambda.eIoni.e-.asc",true);
   
   GPPhysicsTable msc_table;
   ReadPhysicsTable(msc_table,"data/Lambda.msc.e-.asc",true);
   
   printf("Copying GPPhysicsTable data from host to device\n");
   broker->Upload_eBremTable(eBrem_table);
   broker->Upload_eIoniTable(eIoni_table);
   broker->UploadMscTable(msc_table);
   
}

bool CoprocessorBroker::CudaSetup(int nblocks, int nthreads, int maxTrackPerThread)
{
   setup(this);

   
   fNblocks = nblocks;
   fNthreads = nthreads;
   fMaxTrackPerThread = maxTrackPerThread;

   //initialize the stream
   for(unsigned int i=0; i < sizeof(fStream)/sizeof(StreamHelper); ++i) {
      fStream[i].CudaSetup(i,nblocks,nthreads,maxTrackPerThread);
      fStream[i].Push(&fHelpers);
   }
   return true;
}

bool CanHandleTrack(GeantTrack &track)
{
   // Currently we can only handle electron, which we pretend are the only
   // particle to have charge -1.
   
   return -1 == track.charge;
}

unsigned int CoprocessorBroker::StreamHelper::TrackToDevice(int tid,
                           GeantTrack **host_track, int *trackin,
                           unsigned int startIdx, unsigned int basketSize)
{
   if (fThreadId == -1) {
      fThreadId = tid;
      fHostTracks = host_track;
   } else if (fThreadId != tid || fHostTracks != host_track) {
      Fatal("TrackToDevice","The stream %p is already assigned to thread %d when requested to do work for thread %d",
            this, fThreadId, tid);
      return 0;
   }

   WorkloadManager *mgr = WorkloadManager::Instance();
   unsigned int start = fNStaged;
   unsigned int count = 0;
   for(unsigned int hostIdx = startIdx;
       fNStaged < fChunkSize && hostIdx < basketSize;
       ++hostIdx ) {

      ++count;

      TGeoBranchArray *path = host_track[trackin[hostIdx]]->path;
//      if (path->GetLevel()>1) {
//         fprintf(stderr,"DEBUG: for %d level is %d\n",trackin[hostIdx],path->GetLevel());
//      }
      //fprintf(stderr,"DEBUG8: See track #%d\n",trackin[hostIdx]);
      TGeoVolume *logVol = path->GetCurrentNode()->GetVolume();
      int volumeIndex = ((GeantVolumeBasket*)logVol->GetField())->GetNumber();

      if (!CanHandleTrack(* host_track[trackin[hostIdx]])) {
         
         fTrackCollection->AddTrack(trackin[hostIdx],
                                    mgr->GetBasketArray()[volumeIndex]);
         trackin[hostIdx] = -1; // forget the track since we already passed it along.
         gPropagator->fTransportOngoing = kTRUE;
         continue;
      }
      fTrackId[fNStaged]  = trackin[hostIdx];
      fTrack[fNStaged].x  = host_track[trackin[hostIdx]]->xpos;
      fTrack[fNStaged].y  = host_track[trackin[hostIdx]]->ypos;
      fTrack[fNStaged].z  = host_track[trackin[hostIdx]]->zpos;
      fTrack[fNStaged].px = host_track[trackin[hostIdx]]->px;
      fTrack[fNStaged].py = host_track[trackin[hostIdx]]->py;
      fTrack[fNStaged].pz = host_track[trackin[hostIdx]]->pz;
      fTrack[fNStaged].q  = host_track[trackin[hostIdx]]->charge;
      fTrack[fNStaged].s  = host_track[trackin[hostIdx]]->step;
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

      ++fNStaged;
   }
   mgr->CollectorQueue()->push(fTrackCollection);
   fTrackCollection = GetNextTrackCollection();

   // count = min(fChunkSize,basketSize-startIdx);
   int ntrack = fNStaged - start;
   cudaMemcpyAsync(fDevTrack+start, fTrack+start, ntrack*sizeof(GXTrack),
                   cudaMemcpyHostToDevice, fStream);
   cudaMemcpyAsync(fDevTrackLogIndex+start, fLogIndex+start, ntrack*sizeof(int),
                   cudaMemcpyHostToDevice, fStream);
   cudaMemcpyAsync(fDevTrackPhysIndex+start, fPhysIndex+start, ntrack*sizeof(int),
                   cudaMemcpyHostToDevice, fStream);
   
   return count;
}


unsigned int CoprocessorBroker::StreamHelper::TrackToHost()
{
   WorkloadManager *mgr = WorkloadManager::Instance();
   for(unsigned int devIdx = 0;
       devIdx < fNStaged;
       ++devIdx) {
      if ( fTrackId[devIdx] < 0 ) {
         // fprintf(stderr,"DEBUG: missig a track %d %d\n",devIdx,fTrackId[devIdx]);
         continue;
      }
      fHostTracks[fTrackId[devIdx]]->xpos = fTrack[devIdx].x;
      fHostTracks[fTrackId[devIdx]]->ypos = fTrack[devIdx].y;
      fHostTracks[fTrackId[devIdx]]->zpos = fTrack[devIdx].z;
      fHostTracks[fTrackId[devIdx]]->px = fTrack[devIdx].px;
      fHostTracks[fTrackId[devIdx]]->py = fTrack[devIdx].py;
      fHostTracks[fTrackId[devIdx]]->pz = fTrack[devIdx].pz;
      fHostTracks[fTrackId[devIdx]]->charge = fTrack[devIdx].q;
      fHostTracks[fTrackId[devIdx]]->step = fTrack[devIdx].s;
      
      // NOTE: need to update the path ....

      if (fLogIndex[devIdx] == -1) {
         // The track/particle is done.
         // fprintf(stderr,"DEBUG: stopping %d in basket %d %d:%d - isalive %d\n", fTrackId[devIdx],fLogIndex[devIdx],
         //         fHostTracks[fTrackId[devIdx]]->event,fTrackId[devIdx],fHostTracks[fTrackId[devIdx]]->IsAlive());
         gPropagator->StopTrack(gPropagator->fTracks[fTrackId[devIdx]]);
      } else if (1) {
         // Let's reschedule it.
         // fprintf(stderr,"DEBUG: rescheduling %d in basket %d\n",fTrackId[devIdx],fLogIndex[devIdx]);
         fTrackCollection->AddTrack(fTrackId[devIdx],
                                                  mgr->GetBasketArray()[fLogIndex[devIdx]]);
         gPropagator->fTransportOngoing = kTRUE;
      }
   }
   mgr->CollectorQueue()->push(fTrackCollection);
   fTrackCollection = GetNextTrackCollection();
   fThreadId = -1;
   fHostTracks = 0;
   return fNStaged;
}

// typedef void(CUDART_CB * cudaStreamCallback_t)(cudaStream_t stream, cudaError_t status, void
//                                                *userData)

void StreamReset(cudaStream_t /* stream */, cudaError_t status, void *userData)
{
   CoprocessorBroker::StreamHelper *helper = (CoprocessorBroker::StreamHelper*)userData;
   helper->Reset();
}

void TrackToHost(cudaStream_t /* stream */, cudaError_t status, void *userData)
{
   CoprocessorBroker::StreamHelper *helper = (CoprocessorBroker::StreamHelper*)userData;
   helper->TrackToHost();
}


CoprocessorBroker::Stream CoprocessorBroker::GetNextStream()
{
   // Return the current stream (one that we can still add input to)
   // or return new one.

   if (!fCurrentHelper) {
      fCurrentHelper = dynamic_cast<StreamHelper*>(fHelpers.wait_and_pop());
      if (!fCurrentHelper) {
         // nothing we can do at the moment
         return 0;
      }
   }
   return fCurrentHelper;
}

CoprocessorBroker::Stream CoprocessorBroker::launchTask(bool wait /* = false */)
{
   if (fCurrentHelper == 0 || fCurrentHelper->fNStaged == 0) {
      // Nothing to do ...
      return 0;
   }
   StreamHelper *stream = fCurrentHelper;
   fCurrentHelper = 0;
   
   if (stream) {
      Printf("(%d - GPU) == Starting kernel on stream %d with %d tracks\n",
             stream->fThreadId, stream->fStreamId, stream->fNStaged );
      fTotalWork += stream->fNStaged;
      tracking_gpu(stream->fdRandStates,(GPGeomManager*)fdGeometry,fdFieldMap,
                   stream->fDevTrack,
                   stream->fDevTrackLogIndex,
                   stream->fDevTrackPhysIndex,
                   fd_eBremTable, fd_eIoniTable, fd_mscTable,
                   stream->fNStaged,
                   fNblocks,fNthreads,*stream);
      
      cudaMemcpyAsync(stream->fTrack, stream->fDevTrack, stream->fNStaged*sizeof(GXTrack),
                      cudaMemcpyDeviceToHost, *stream);
      cudaMemcpyAsync(stream->fPhysIndex, stream->fDevTrackPhysIndex, stream->fNStaged*sizeof(int),
                      cudaMemcpyDeviceToHost, *stream);
      cudaMemcpyAsync(stream->fLogIndex, stream->fDevTrackLogIndex, stream->fNStaged*sizeof(int),
                      cudaMemcpyDeviceToHost, *stream);
      
      cudaStreamAddCallback(stream->fStream, TrackToHost, stream, 0 );
      cudaStreamAddCallback(stream->fStream, StreamReset, stream, 0 );
      
      if (wait) {
         cudaStreamSynchronize(*stream);
//         stream->ResetNStaged();
//         fHelpers.push(stream);
      }
   }
   return stream;
}

void CoprocessorBroker::runTask(int threadid, int nTracks, int volumeIndex, GeantTrack **tracks, int *trackin)
{
   bool force = false;
   
   
   //<<<---------------------------------------------------------->>>
   if(fKernelType ==1) {
//         fprintf(stderr,"round = %d count1 = %d fNchunk=%d nTracks=%d\n",
//                 round,count1,fNchunk,nTracks);
//         cudaMemcpyAsync(track_d0, track_h+i, count1*sizeof(GXTrack),
//                         cudaMemcpyHostToDevice, stream0);
//
      int trackLeft  = nTracks;
      int trackStart = 0;
      while (trackLeft) {
         StreamHelper *stream = GetNextStream();
         if (!stream) return;

         cudaEvent_t start, stop;
         cudaEventCreate (&start);
         cudaEventCreate (&stop);
         
         //start time for kernels
         cudaEventRecord (start,0);

         unsigned int before = stream->fNStaged;
         unsigned int count = stream->TrackToDevice(threadid,
                                                    tracks, trackin,
                                                    trackStart, nTracks);
         unsigned int rejected = nTracks-trackStart - (stream->fNStaged-before);
         if (stream->fNStaged < stream->fChunkSize
             && (before != stream->fNStaged) // if we did not make any progress, assume there is no 'interesting' track left and schedule the kernel.
             && !force) {
            // We do not have enough tracks
            Printf("(%d - GPU) ================= Stream %d Tracks: %d seen %d skipped %d accumulated", threadid, stream->fStreamId, count, rejected, stream->fNStaged);
            
            // Stop time for kernel
            cudaEventRecord (stop,0);
            cudaEventSynchronize (stop);
            float elapsedTimeGPU;
            cudaEventElapsedTime (&elapsedTimeGPU,start,stop);

            return;
         }

         if (stream->fNStaged) {
            Printf("(%d - GPU) ================= Stream %d Tracks: %d seen %d skipped %d accumulated", threadid, stream->fStreamId, count, rejected, stream->fNStaged);
            launchTask();
         }
         trackLeft  -= count;
         trackStart += count;

         cudaEventRecord (stop,0);
         cudaEventSynchronize (stop);
         float elapsedTimeGPU;
         cudaEventElapsedTime (&elapsedTimeGPU,start,stop);
      }
   }
   else if(fKernelType ==2) {
      
      fprintf(stderr,"We do not yet support the multiple kernel\n");
#if 0
      //photon kernel
      for(int i = 0 ; i < nPhotons ; i += 2*fNchunkG) {
         int count1 = min(fNchunkG,nPhotons-i) ;
         cudaMemcpyAsync(photon_d0, photon_h+i, count1*sizeof(GXTrack),
                         cudaMemcpyHostToDevice, stream0);
         
         tracking_photon_gpu(fdRandStates0,(GPVPhysicalVolume*)fdGeometry,fdFieldMap,photon_d0,
                             fd_eBremTable, fd_eIoniTable, fd_mscTable,
                             count1,
                             theNBlocks,theNThreads, stream0);
         
         cudaMemcpyAsync(photon_h+i, photon_d0, count1*sizeof(GXTrack),
                         cudaMemcpyDeviceToHost, stream0);
         
         cudaMemcpyAsync(photon_d1, photon_h+i+fNchunkG, (nPhotons/nDiv)*sizeof(GXTrack),
                         cudaMemcpyHostToDevice, stream1);
         
         tracking_photon_gpu(fdRandStates1,(GPVPhysicalVolume*)fdGeometry,fdFieldMap,photon_d1,
                             fd_eBremTable, fd_eIoniTable, fd_mscTable,
                             (nPhotons/nDiv), theNBlocks,theNThreads, stream1);
         
         cudaMemcpyAsync(photon_h+i+fNchunkG, photon_d1, (nPhotons/nDiv)*sizeof(GXTrack),
                         cudaMemcpyDeviceToHost, stream1);
      }
      
      //electon kernel
      for(int i = 0 ; i < nElectrons ; i += 2*fNchunkE) {
         int count1 = min(fNchunkE,nElectrons-i) ;
         
         cudaMemcpyAsync(electron_d0, electron_h+i, count1*sizeof(GXTrack),
                         cudaMemcpyHostToDevice, stream0);
         
         tracking_electron_gpu(fdRandStates0,(GPVPhysicalVolume*)fdGeometry,fdFieldMap,electron_d0,
                               fd_eBremTable, fd_eIoniTable, fd_mscTable,
                               count1,
                               theNBlocks,theNThreads, stream0);
         
         cudaMemcpyAsync(electron_h+i, electron_d0, count1*sizeof(GXTrack),
                         cudaMemcpyDeviceToHost, stream0);
         
         int count2 = min(fNchunkE, nElectrons - (i+count1));
         cudaMemcpyAsync(electron_d1, electron_h+i+count1, count2*sizeof(GXTrack),
                         cudaMemcpyHostToDevice, stream1);
         
         tracking_electron_gpu(fdRandStates1,(GPVPhysicalVolume*)fdGeometry,fdFieldMap,electron_d1,
                               fd_eBremTable, fd_eIoniTable, fd_mscTable,
                               count2,
                               theNBlocks,theNThreads, stream1);
         
         cudaMemcpyAsync(electron_h+i+count1, electron_d1, count2*sizeof(GXTrack),
                         cudaMemcpyDeviceToHost, stream1);
      }
#endif
   }
   
}

void CoprocessorBroker::waitForTasks()
{
   // Make sure all the tasks are finished.
   
   for(unsigned int i=0; i < sizeof(fStream)/sizeof(StreamHelper); ++i) {
      if (fStream[i].fNStaged) {
         cudaStreamSynchronize(fStream[i]);
      }
   }
}
