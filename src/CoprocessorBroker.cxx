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


CoprocessorBroker::StreamHelper::StreamHelper() : fTrack(0),fPhysIndex(0),fLogIndex(0),fNStaged(0)
{
   // Default constructor.
}

CoprocessorBroker::StreamHelper::~StreamHelper() {
   // Destructor.
   
   delete [] fTrack;
   delete [] fPhysIndex;
   delete [] fLogIndex;
   cudaStreamDestroy(fStream);
}


bool CoprocessorBroker::StreamHelper::CudaSetup(int nblocks, int nthreads, int maxTrackPerThread)
{
   cudaStreamCreate(&fStream);
   
   //prepare random engines on the device
   fdRandStates.Alloc( nblocks*nthreads );
   curand_setup_gpu(fdRandStates, time(NULL), nblocks, nthreads);
   
   unsigned int maxTrackPerKernel = nblocks*nthreads*maxTrackPerThread;
   fDevTrack.Alloc(maxTrackPerKernel);
   fDevTrackPhysIndex.Alloc(maxTrackPerKernel);
   fDevTrackLogIndex.Alloc(maxTrackPerKernel);
   fTrack = new GXTrack[maxTrackPerKernel];
   fPhysIndex = new int[maxTrackPerKernel];
   fLogIndex = new int[maxTrackPerKernel];
   
   return true;
}

CoprocessorBroker::CoprocessorBroker() : fdGeometry(0)
   ,fNblocks(0),fNthreads(0),fNchunk(0)
   ,fMaxTrackPerThread(1)
   ,fKernelType(1)
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
   fStream0.CudaSetup(nblocks,nthreads,maxTrackPerThread);
   fStream1.CudaSetup(nblocks,nthreads,maxTrackPerThread);

   return true;
}

void CoprocessorBroker::prepateDataArray(long nTracks)
{
   
   if(fKernelType ==1) {
      // Nchunk = ( nTracks%nDiv != 0 ) ? (nTracks/nDiv + 1 ) : (nTracks/nDiv);
      int nthreads_total = fNblocks * fNthreads;
      if (nTracks < nthreads_total*fMaxTrackPerThread) {
         fNchunk = nTracks;
         fprintf(stderr,"Warning: too few tracks resulting in only one step with low chunk size %d (< %d)\n",fNchunk,nthreads_total);
      } else {
         fNchunk = nthreads_total*fMaxTrackPerThread;
         // if ( (nTracks%fNchunk)==0 )
         //   nstep = nTracks / fNchunk;
         // else
         //   nstep = nTracks / fNchunk + 1;
      }

//      cudaMalloc((void**)&track_d0, Nchunk*sizeof(GXTrack));
//      cudaMalloc((void**)&track_d1, Nchunk*sizeof(GXTrack));
      //      cudaMalloc((void**)&track_d1, Nchunk*sizeof(GXTrack));
   }
   else if(fKernelType ==2) {
      fprintf(stderr,"We do not yet support the multiple kernel\n");
#if 0
      
      int nthreads_total = theNBlocks * theNThreads;
      
      NchunkG = (nPhotons/nDiv);
      if (NchunkG >= nthreads_total) {
         NchunkG -= ( (nPhotons/nDiv) % nthreads_total); // align on multiple of the grid size
      } else {
         fprintf(stderr,"Warning: too many division for photons (%d) resulting in low chunk size %d (< %d)\n",nDiv,NchunkG,nthreads_total);
      }
      
      NchunkE = (nElectrons/nDiv);
      if (NchunkE >= nthreads_total) {
         NchunkE -= ( (nElectrons/nDiv) % nthreads_total); // align on multiple of the grid size
      } else {
         fprintf(stderr,"Warning: too many division for electrons (%d) resulting in low chunk size %d (< %d)\n",nDiv,NchunkE,nthreads_total);
      }
#endif
//      cudaMalloc((void**)&photon_d0, NchunkG*sizeof(GXTrack));
//      cudaMalloc((void**)&photon_d1, NchunkG*sizeof(GXTrack));
//      
//      cudaMalloc((void**)&electron_d0, NchunkE*sizeof(GXTrack));
//      cudaMalloc((void**)&electron_d1, NchunkE*sizeof(GXTrack));
   }
   
}

bool CanHandleTrack(GeantTrack &track)
{
   // Currently we can only handle electron, which we pretend are the only
   // particle to have charge -1.
   
   return -1 == track.charge;
}

unsigned int CoprocessorBroker::StreamHelper::TrackToDevice(int tid,
                           GeantTrack **host_track, int *trackin,
                           unsigned int startIdx, unsigned int basketSize,
                           unsigned int chunkSize)
{
   unsigned int count = 0;
   for(unsigned int hostIdx = startIdx;
       fNStaged < chunkSize && hostIdx < basketSize;
       ++hostIdx ) {

      ++count;

      TGeoBranchArray *path = host_track[trackin[hostIdx]]->path;
//      if (path->GetLevel()>1) {
//         fprintf(stderr,"DEBUG: for %d level is %d\n",trackin[hostIdx],path->GetLevel());
//      }
      //fprintf(stderr,"DEBUG8: See track #%d\n",trackin[hostIdx]);
      TGeoVolume *logVol = path->GetCurrentNode()->GetVolume();
      int volumeIndex = ((GeantVolumeBasket*)logVol->GetField())->GetNumber();

      if (1 && !CanHandleTrack(* host_track[trackin[hostIdx]])) {
         
         gPropagator->fCollections[tid]->AddTrack(trackin[hostIdx],
                                                  WorkloadManager::Instance()->GetBasketArray()[volumeIndex]);
         trackin[hostIdx] = -1; // forget the track since we already passed it along.
         gPropagator->fTransportOngoing = kTRUE;
         continue;
      }
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
   // count = min(chunkSize,basketSize-startIdx);
   cudaMemcpyAsync(fDevTrack, fTrack, count*sizeof(GXTrack),
                   cudaMemcpyHostToDevice, fStream);
   cudaMemcpyAsync(fDevTrackLogIndex, fLogIndex, count*sizeof(int),
                   cudaMemcpyHostToDevice, fStream);
   cudaMemcpyAsync(fDevTrackPhysIndex, fPhysIndex, count*sizeof(int),
                   cudaMemcpyHostToDevice, fStream);
   
   return count;
}


unsigned int CoprocessorBroker::StreamHelper::TrackToHost(// int tid,
                                                            GeantTrack **host_track, int *trackin,
                                                            unsigned int startIdx)
{
   for(unsigned int devIdx = 0, hostIdx = startIdx;
       devIdx < fNStaged;
       ++devIdx, ++hostIdx ) {
      if ( trackin[hostIdx] < 0 ) {
         --devIdx;
         continue;
      }
      host_track[trackin[hostIdx]]->xpos = fTrack[devIdx].x;
      host_track[trackin[hostIdx]]->ypos = fTrack[devIdx].y;
      host_track[trackin[hostIdx]]->zpos = fTrack[devIdx].z;
      host_track[trackin[hostIdx]]->px = fTrack[devIdx].px;
      host_track[trackin[hostIdx]]->py = fTrack[devIdx].py;
      host_track[trackin[hostIdx]]->pz = fTrack[devIdx].pz;
      host_track[trackin[hostIdx]]->charge = fTrack[devIdx].q;
      host_track[trackin[hostIdx]]->step = fTrack[devIdx].s;
      
      // NOTE: need to update the path ....
   }
   return fNStaged;
}

struct TrackToHostData {
   GeantTrack **host_track;
   Int_t       *trackin;
   UInt_t       startIdx;
   CoprocessorBroker::StreamHelper *fStream;
   
   TrackToHostData(GeantTrack **i_host_track, Int_t *i_trackin, UInt_t i_startIdx,
                   CoprocessorBroker::StreamHelper *stream) :
   host_track(i_host_track), trackin(i_trackin), startIdx(i_startIdx),
   fStream(stream)
   {}

   unsigned int TrackToHost()
   {
      return fStream->TrackToHost(host_track,trackin,startIdx);
   }
   
};


// typedef void(CUDART_CB * cudaStreamCallback_t)(cudaStream_t stream, cudaError_t status, void
//                                                *userData)

void ResetNStaged(cudaStream_t /* stream */, cudaError_t status, void *userData)
{
   CoprocessorBroker::StreamHelper *helper = (CoprocessorBroker::StreamHelper*)userData;
   helper->ResetNStaged();
}

void TrackToHost(cudaStream_t /* stream */, cudaError_t status, void *userData)
{
   TrackToHostData *data = (TrackToHostData*)userData;
   data->TrackToHost();
}


void CoprocessorBroker::runTask(int threadid, int nTracks, int volumeIndex, GeantTrack **tracks, int *trackin)
{
   
   cudaEvent_t start, stop;
   cudaEventCreate (&start);
   cudaEventCreate (&stop);
 
   //start time for kernels
   cudaEventRecord (start,0);
   
   std::vector<TrackToHostData> data;
   data.reserve( (nTracks / fNchunk) + 1 ); // Important, we can not afford a resize.

   //<<<---------------------------------------------------------->>>
   if(fKernelType ==1) {
      for(int i = 0, round = 0 ; i < nTracks ; i += 2*fNchunk, ++round) {
         
         int count1 = std::min(fNchunk,nTracks-i);
//         fprintf(stderr,"round = %d count1 = %d fNchunk=%d nTracks=%d\n",
//                 round,count1,fNchunk,nTracks);
//         cudaMemcpyAsync(track_d0, track_h+i, count1*sizeof(GXTrack),
//                         cudaMemcpyHostToDevice, stream0);
//
         count1 = fStream0.TrackToDevice(threadid,
                                         tracks,trackin,
                                         i,nTracks,
                                         fNchunk);
         printf("Running kernel 1 with %d tracks\n",count1);
         cudaStreamAddCallback(fStream0.fStream, ResetNStaged, &fStream0, 0 );
         
         tracking_gpu(fStream0.fdRandStates,(GPGeomManager*)fdGeometry,fdFieldMap,
                      fStream0.fDevTrack,
                      fStream0.fDevTrackLogIndex,
                      fStream0.fDevTrackPhysIndex,
                      fd_eBremTable, fd_eIoniTable, fd_mscTable,
                      fStream0.fNStaged,
                      fNblocks,fNthreads,fStream0);
         
         cudaMemcpyAsync(fStream0.fTrack, fStream0.fDevTrack, fStream0.fNStaged*sizeof(GXTrack),
                         cudaMemcpyDeviceToHost, fStream0);
         cudaMemcpyAsync(fStream0.fPhysIndex, fStream0.fDevTrackPhysIndex, fStream0.fNStaged*sizeof(int),
                         cudaMemcpyDeviceToHost, fStream0);
         cudaMemcpyAsync(fStream0.fLogIndex, fStream0.fDevTrackLogIndex, fStream0.fNStaged*sizeof(int),
                         cudaMemcpyDeviceToHost, fStream0);
         
         data.push_back(TrackToHostData(tracks, trackin, i, &fStream0));
         cudaStreamAddCallback(fStream0.fStream, TrackToHost, &(data.back()), 0 );

         int count2 = std::min( fNchunk, nTracks - (i+count1) );
         if (count2) {
//            cudaMemcpyAsync(track_d1, track_h+i+count1,count2*sizeof(GXTrack),
//                            cudaMemcpyHostToDevice, stream1);
            count2 = fStream1.TrackToDevice(threadid, tracks,trackin,i+count1,nTracks,
                                            fNchunk);
           
            printf("Running kernel 1 with %d tracks\n",count2);
            cudaStreamAddCallback(fStream1.fStream, ResetNStaged, &fStream1, 0 );
            tracking_gpu(fStream1.fdRandStates,(GPGeomManager*)fdGeometry,fdFieldMap,
                         fStream1.fDevTrack,
                         fStream0.fDevTrackLogIndex,
                         fStream1.fDevTrackPhysIndex,
                         fd_eBremTable, fd_eIoniTable, fd_mscTable,
                         fStream1.fNStaged,
                         fNblocks,fNthreads,fStream1);
            
            cudaMemcpyAsync(fStream1.fTrack, fStream1.fDevTrack, fStream1.fNStaged*sizeof(GXTrack),
                            cudaMemcpyDeviceToHost, fStream1.fStream);
            cudaMemcpyAsync(fStream1.fPhysIndex, fStream1.fDevTrackPhysIndex, fStream1.fNStaged*sizeof(int),
                            cudaMemcpyDeviceToHost, fStream1);
            cudaMemcpyAsync(fStream1.fLogIndex, fStream1.fDevTrackLogIndex, fStream1.fNStaged*sizeof(int),
                            cudaMemcpyDeviceToHost, fStream1);
            //{userdata} = { host_track, trackin, i+count1, scratch_track, count2 }
            data.push_back(TrackToHostData(tracks, trackin, i+count1, &fStream1));
            cudaStreamAddCallback(fStream1, TrackToHost, &(data.back()), 0 );
         }
         // fprintf(stderr,"round = %d count1 = %d count2 = %d\n",round,count1,count2);
         // random_testing_gpu( devStates, results,theNBlocks,theNThreads,stream0);
         
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
   
   cudaStreamSynchronize(fStream0);
   cudaStreamSynchronize(fStream1);
   
   //stop time for kernel
   cudaEventRecord (stop,0);
   cudaEventSynchronize (stop);
   float elapsedTimeGPU;
   cudaEventElapsedTime (&elapsedTimeGPU,start,stop);
      
}


void retrieveData()
{
   
   
   
}
