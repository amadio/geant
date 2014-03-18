#include <iostream>
#include <fstream>

// #include <cuda.h>
// #include <curand.h>
// #include <cuda_runtime.h>
#include <stdio.h>
#include <math.h>

#include "tracking_kernel.h"
//#include "random_kernel.h"

#include "GPConstants.h"
//#include "GPFieldMap.h"
//#include "GPFieldMapData.h"
#include "cmsExpMagneticField.hh"

#include "G4Track.hh"
#include "CLHEP/Vector/ThreeVector.h"

//Geometry
#include "GPVGeometry.h"
#include "GPUserGeometry.h"
#include "GPSimpleEcal.h"
//#include "GPSimpleCMS.h"
#include "G4VPhysicalVolume.hh"

//EMPhysics
//#include "GPConstants.h"
#include "GXEMPhysicsUtils.h"

#include "G4ParticleTable.hh"
#include "G4Electron.hh"
#include "G4OpticalPhoton.hh"
#include "stopwatch.h"

#include "G4ProductionCutsTable.hh"
#include "G4Region.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"

using namespace std;

int main (int argc, char* argv[]) 
{
  //argument
  int kernelType = 1; //1. general traportation kernel
  //2. electron/photon kernels
  int flagTest = 0; //1 test mode
  // int ndivarg = -1;

  if(argc >= 2) kernelType = atoi(argv[1]);
  if(argc >= 3) flagTest = atoi(argv[2]);
  // if(argc >= 4) ndivarg = atoi(argv[3]);

  if(kernelType >2 ) {
    std::cout << "Wrong kernelType: Bailing out ... "<< std::endl;
    return 0;
  }
  
  //check whether device overlap is available
  // int whichDevice;
  // cudaGetDevice(&whichDevice);
  
  // if(whichDevice > 0) {
  //   cudaDeviceReset();
    
  //   cudaDeviceProp prop;
  //   cudaGetDeviceProperties(&prop,whichDevice);
  //   if(!prop.deviceOverlap) {
  //     printf("Device overlap is not available\n");
  //     return 0;
  //   }
  // }
  
  //default number of blocks and threads
  
  // int theNBlocks  =  32;
  // int theNThreads = 128;
  //int theNBlocks  = 128;
  //int theNThreads = 448;
  
  //reconfigure 
  // if(flagTest>0) {
  //   theNBlocks = 2;
  //   theNThreads = 10;
  // }
  
  // char* cudaNBlocks = getenv("GP_CUDA_NBLOCKS");
  // char* cudaNThreads = getenv("GP_CUDA_NTHREADS");
  
  // if(cudaNBlocks) theNBlocks = atoi(cudaNBlocks);
  // if(cudaNThreads) theNThreads = atoi(cudaNThreads);
  
  // std::cout << "... tracking_kernel<<<" << theNBlocks << "," 
  // << theNThreads <<">>> (...) ..." << std::endl;
  
  //1. Construct geometry
  int nphi = 4;
  int nz   = 3;
  double density = 8.28;
  
  // From G4RunManagerKernel.cc
  GPVGeometry *geom = new GPSimpleEcal(nphi,nz,density);
  //  GPVGeometry *geom = new GPSimpleCMS();
  
  geom->create();

  G4Region *defaultRegion = new G4Region("DefaultRegionForTheWorld"); // deleted by store
  defaultRegion->UsedInMassGeometry(true);
  defaultRegion->SetProductionCuts(
    G4ProductionCutsTable::GetProductionCutsTable()->GetDefaultProductionCuts());
  geom->getWorldLogicalVolume()->SetRegion(defaultRegion);
  defaultRegion->SetWorld(geom->getWorldVolume());
  defaultRegion->AddRootLogicalVolume(geom->getWorldLogicalVolume());
  defaultRegion->UpdateMaterialList();
  
  G4ProductionCutsTable::GetProductionCutsTable()->UpdateCoupleTable(geom->getWorldVolume());
  G4ProductionCutsTable::GetProductionCutsTable()->RetrieveCutsTable("./data",true); 

  //Navigator
   G4Navigator *aNavigator = new G4Navigator();
  aNavigator->SetWorldVolume(geom->getWorldVolume());
  //G4Navigator_Constructor(&aNavigator);
  //G4Navigator_SetWorldVolume(&aNavigator,world);
   
  G4TransportationManager::GetTransportationManager()->SetNavigatorForTracking(aNavigator);

  // GPVGeometry::byte *geom_h = (GPVGeometry::byte *) malloc (geom->size()) ;
  // geom->relocate( geom_h );
  // memcpy(geom_h,geom->getBuffer(),geom->size());  
  
  // GPVGeometry::byte *geom_d;
  // // cudaMalloc( (void**)&geom_d, geom->size() );
  // geom->relocate( geom_d );
  // // cudaMemcpy(geom_d, geom->getBuffer(), geom->size(), cudaMemcpyHostToDevice);
  
  //2. Read magnetic field map

  // GPFieldMap *bmap_d;
  // GPFieldMap *bmap_h;
  // GPFieldMap** fieldMap;
  
  // fieldMap = (GPFieldMap **) malloc (nbinZ*sizeof(GPFieldMap *));
  // for (int j = 0 ; j < nbinZ ; j++) {
  //   fieldMap[j] = (GPFieldMap *) malloc (nbinR*sizeof(GPFieldMap));
  // } 
  
  const char* fieldMapFile = getenv("GP_BFIELD_MAP");
  fieldMapFile = (fieldMapFile) ? fieldMapFile : "cmsExp.mag.3_8T";
  
  cmsExpMagneticField *bmap = new cmsExpMagneticField();
  bmap->ReadFieldMap( fieldMapFile );
  bmap->SetFieldType("volumebase");

  // std::ifstream ifile(fieldMapFile, ios::in | ios::binary | ios::ate);
  
  // if (ifile.is_open()) {
    
  //   //field map structure
  //   GPFieldMapData fd;
    
  //   ifstream::pos_type fsize = ifile.tellg();
  //   size_t dsize = sizeof(GPFieldMapData);    
    
  //   long int ngrid = fsize/dsize;
  //   ifile.seekg (0, ios::beg);
    
  //   std::cout << "... transportation ... Loading magnetic field map: " 
  //   << fieldMapFile << std::endl;
    
  //   for(int i = 0 ; i < ngrid ; i++) {
  //     ifile.read((char *)&fd, sizeof(GPFieldMapData));
      
  //     //check validity of input data
  //     if(abs(fd.iz) > noffZ || fd.ir > nbinR) {
  //       std::cout << " Field Map Array Out of Range" << std::endl;
  //     }
  //     else {
  //       fieldMap[noffZ+fd.iz][fd.ir].Bz = fd.Bz; 
  //       fieldMap[noffZ+fd.iz][fd.ir].Br = fd.Br;
  //     }
  //   }
  //   ifile.close();
  // }
  
  //3. Create magnetic field on the device
  
  // printf("Creating magnetic field map on the GPU device\n");
  // //prepare fieldMap array: fieldMap[nbinZ][nbinR];
  
  // bmap_h = (GPFieldMap *) malloc (nbinZ*nbinR*sizeof(GPFieldMap));
  
  // for (int i = 0 ; i < nbinZ ; i++) {
  //   for (int j = 0 ; j < nbinR ; j++) {
  //     bmap_h[i+j*nbinZ].Bz = fieldMap[i][j].Bz;
  //     bmap_h[i+j*nbinZ].Br = fieldMap[i][j].Br;
  //   }
  // }
  
  // printf("Copying data from host to device\n");
  
  // cudaMalloc((void**)&bmap_d, nbinZ*nbinR*sizeof(GPFieldMap)) ;
  // cudaMemcpy(bmap_d,bmap_h,nbinZ*nbinR*sizeof(GPFieldMap),
  //            cudaMemcpyHostToDevice);
  
  // 4. Prepare EM physics tables
  bool useSpline = true;

  // We should really test for the presence of the files!

  G4PhysicsTable *eBrem_table = new G4PhysicsTable();
  readTable(eBrem_table,"data/Lambda.eBrem.e-.asc");
  int nv = eBrem_table->length(); // nPhysicsVector;
  for(int j=0; j < nv; j++){
     // eBrem_table.physicsVectors[j].SetSpline(useSpline);
     (*eBrem_table)(j)->SetSpline(useSpline);
  }
  //  GPPhysicsTable_Print(&eBrem_table);
  
  G4PhysicsTable *eIoni_table = new G4PhysicsTable();
  readTable(eIoni_table,"data/Lambda.eIoni.e-.asc");
  nv = eIoni_table->length(); // .nPhysicsVector;
  for(int j=0; j < nv; j++){
     (*eIoni_table)(j)->SetSpline(useSpline);
  }
  //  GPPhysicsTable_Print(&eIoni_table);
  
  G4PhysicsTable * msc_table = 0;
//  readTable(&msc_table,"data/Lambda.msc.e-.asc");
//  nv = msc_table.length(); // .nPhysicsVector;
//  for(int j=0; j < nv; j++){
//     msc_table(j)->SetSpline(useSpline);
//  }
  //  GPPhysicsTable_Print(&msc_table);
  
  printf("Copying GPPhysicsTable data from host to device\n");
  
  // GPPhysicsTable* eBrem_table_d;
  // GPPhysicsTable* eIoni_table_d;
  // GPPhysicsTable* msc_table_d;
  
  // cudaMalloc((void**)&eBrem_table_d, sizeof(GPPhysicsTable));
  // cudaMemcpy(eBrem_table_d, &eBrem_table, sizeof(GPPhysicsTable),cudaMemcpyHostToDevice);
  
  // cudaMalloc((void**)&eIoni_table_d, sizeof(GPPhysicsTable));
  // cudaMemcpy(eIoni_table_d, &eIoni_table, sizeof(GPPhysicsTable),cudaMemcpyHostToDevice);
  
  // cudaMalloc((void**)&msc_table_d, sizeof(GPPhysicsTable));
  // cudaMemcpy(msc_table_d, &msc_table, sizeof(GPPhysicsTable),cudaMemcpyHostToDevice);
  
  //5. Get the list of tracks from the stack
  printf("Preparing a bundle tracks at host\n");
  
  //event-loop
  int nevent = 10;
  
  //initialize the stream
  // cudaStream_t stream0;
  // cudaStream_t stream1;
  // cudaStreamCreate(&stream0);
  // cudaStreamCreate(&stream1);
  
  
  size_t nTracks = 4096*24;
  //int nTracks = 100000;
  // G4Track *track_h;
  G4Track **track_c = new G4Track*[nTracks];
  
  // cudaHostAlloc((void**) &track_h, nTracks*sizeof(GXTrack),
  //               cudaHostAllocDefault);
  // cudaHostAlloc((void**) &track_c, nTracks*sizeof(GXTrack),
                // cudaHostAllocDefault);
  
  // G4Track *electron_h;
  G4Track **electron_c = new G4Track*[nTracks];
  
  // G4Track *photon_h;
  G4Track **photon_c = new G4Track*[nTracks];
  
  if (kernelType == 2) {
    // int nPhotons = nTracks;
    // int nElectrons = nTracks;
    
    // cudaHostAlloc((void**) &photon_h, nPhotons*sizeof(GXTrack),
    //               cudaHostAllocDefault);
    // cudaHostAlloc((void**) &photon_c, nPhotons*sizeof(GXTrack),
    //               cudaHostAllocDefault);
    
    // cudaHostAlloc((void**) &electron_h, nElectrons*sizeof(GXTrack),
    //               cudaHostAllocDefault);
    // cudaHostAlloc((void**) &electron_c, nElectrons*sizeof(GXTrack),
    //               cudaHostAllocDefault);   
  }
  // double *results = 0;
  // cudaMalloc((void**)&results, sizeof(double)*nTracks);
  
  for (int i = 0 ; i < nevent ; i++) {
    
    // int nDiv = ndivarg > 0 ? ndivarg : 24;
    
    // int Nchunk  =       0;
    // int NchunkG  =      0;
    // int NchunkE  =      0;
    
    if(flagTest>0) {
      nTracks = 20;
    }
    
    printf("Event %d> Transportation in kernel for %lu tracks\n",i,nTracks);
    
    //populate tracks with track position, momentum, energy, steplength
    
    int nPhotons = 0;    
    int nElectrons = 0;    
    
    //photon probability: default 20%
    G4double photonFraction = 0.2;
    
    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* photon;
    //photon = particleTable->FindParticle("opticalphoton"); 
    photon = G4OpticalPhoton::OpticalPhotonDefinition();
    G4ParticleDefinition* electron ;
    //electron = particleTable->FindParticle("e-"); 
    electron = G4Electron::ElectronDefinition();
   
    G4Step *steparr = new G4Step[nTracks];

    //randomize position and momentum    
    for(size_t i = 0 ; i < nTracks ; i++){
       //barrel only
       // track_h[i].x      = track_c[i].x     = 300*(2.0*rand()/RAND_MAX-1.0);
       // track_h[i].y      = track_c[i].y     = 300*(2.0*rand()/RAND_MAX-1.0);
       // track_h[i].z      = track_c[i].z     = 300*(2.0*rand()/RAND_MAX-1.0);
       G4ThreeVector position(300*(2.0*rand()/RAND_MAX-1.0),300*(2.0*rand()/RAND_MAX-1.0),300*(2.0*rand()/RAND_MAX-1.0));
       // track_h[i].s      = track_c[i].s     = 1.0+1*(2.0*rand()/RAND_MAX-1.0);
       G4double steplen = 1.0+1*(2.0*rand()/RAND_MAX-1.0); // is 's' the same as step length?
                              
       // track_h[i].px     = track_c[i].px    = 1.+1000*(2.0*rand()/RAND_MAX-1.0);
       // track_h[i].py     = track_c[i].py    = 1.+1000*(2.0*rand()/RAND_MAX-1.0);
       // track_h[i].pz     = track_c[i].pz    = 1.+1000*(2.0*rand()/RAND_MAX-1.0);
       G4ThreeVector momentum( G4ThreeVector(1.+1000*(2.0*rand()/RAND_MAX-1.0),1.+1000*(2.0*rand()/RAND_MAX-1.0),1.+1000*(2.0*rand()/RAND_MAX-1.0) ) );
       G4double mag = momentum.mag();


       // track_c[i].SetMomentumDirection( momentum.unit() );
       // track_c[i].SetKineticEnergy(mag);
       // ((G4DynamicParticle*)track_c[i].GetDynamicParticle())->SetCharge( -1.0 );
       // track_h[i].q      = track_c[i].q     = -1.0;
      
       G4DynamicParticle *particle;
       if( 1.0*rand()/(float)RAND_MAX < photonFraction) {
         // track_h[i].q      = track_c[i].q     = 0.0;
         // ((G4DynamicParticle*)track_c[i].GetDynamicParticle())->SetCharge(0.0);
          particle = new G4DynamicParticle(photon, momentum.unit(), mag);
          nPhotons++;
       }
       else {
          // track_h[i].q      = track_c[i].q     = -1.0;
          particle = new G4DynamicParticle(electron, momentum.unit(), mag);
          nElectrons++;
       }
       track_c[i] = new G4Track(particle, 0.0, position);
       // track_c[i]->SetMomentumDirection( momentum.unit() );
       track_c[i]->SetKineticEnergy( mag );
       track_c[i]->SetStepLength( steplen );
       track_c[i]->SetStep(&(steparr[i]));
    }
    
    if(kernelType ==2) {
      int j = 0;
      int k = 0;
      
      for(size_t i = 0 ; i < nTracks ; i++){
        if(track_c[i]->GetDynamicParticle()->GetCharge() == 0.0) {
          // photon_h[j].x      = photon_c[j].x     = track_h[i].x  ;
          // photon_h[j].y      = photon_c[j].y     = track_h[i].y  ;
          // photon_h[j].z      = photon_c[j].z     = track_h[i].z  ;
          // photon_h[j].s      = photon_c[j].s     = track_h[i].s  ;
          
          // photon_h[j].px     = photon_c[j].px    = track_h[i].px ;
          // photon_h[j].py     = photon_c[j].py    = track_h[i].py ;
          // photon_h[j].pz     = photon_c[j].pz    = track_h[i].pz ;
          // photon_h[j].q      = photon_c[j].q     = track_h[i].q ;
           photon_c[j] = track_c[i]; // .CopyTrackInfo( track_c[i] );
          j++;  
        }
        else {
          // electron_h[k].x      = electron_c[k].x     = track_h[i].x  ;
          // electron_h[k].y      = electron_c[k].y     = track_h[i].y  ;
          // electron_h[k].z      = electron_c[k].z     = track_h[i].z  ;
          // electron_h[k].s      = electron_c[k].s     = track_h[i].s  ;
          
          // electron_h[k].px     = electron_c[k].px    = track_h[i].px ;
          // electron_h[k].py     = electron_c[k].py    = track_h[i].py ;
          // electron_h[k].pz     = electron_c[k].pz    = track_h[i].pz ;
          // electron_h[k].q      = electron_c[k].q     = track_h[i].q ;
           electron_c[k] = track_c[i]; //.CopyTrackInfo( track_c[i] );
          k++;
        }
      }
    }
    
    printf(" (nElectron,nPhoton) = (%d,%d)\n",nElectrons,nPhotons);
    
    // float elapsedTimeAlloc = 0.0;
    // float elapsedTimeGPU = 0.0;
    
    // printf("Copying track data from host to device\n");
    
    //start time memory allocation on GPU
    // cudaEvent_t start, stop;
    // cudaEventCreate (&start);
    // cudaEventCreate (&stop);
    // cudaEventRecord (start,0);
    
    // G4Track *track_d0;
    // G4Track *track_d1;
    
    // G4Track *photon_d0;
    // G4Track *photon_d1;
    
    // G4Track *electron_d0;
    // G4Track *electron_d1;
    
    // if(kernelType ==1) {
    //    // Nchunk = ( nTracks%nDiv != 0 ) ? (nTracks/nDiv + 1 ) : (nTracks/nDiv);
    //    int nthreads_total = theNBlocks * theNThreads;
    //    Nchunk = (nTracks/nDiv);
    //    if (Nchunk >= nthreads_total) {
    //       Nchunk -= ( (nTracks/nDiv) % nthreads_total); // align on multiple of the grid size
    //    } else {
    //       fprintf(stderr,"Warning: too many division (%d) resulting in low chunk size %d (< %d)\n",nDiv,Nchunk,nthreads_total);
    //    }
       
    //    // cudaMalloc((void**)&track_d0, Nchunk*sizeof(GXTrack));
    //    // cudaMalloc((void**)&track_d1, Nchunk*sizeof(GXTrack));
    //    //      cudaMalloc((void**)&track_d1, Nchunk*sizeof(GXTrack));
    // } 
    // else if(kernelType ==2) {
       
    //    int nthreads_total = theNBlocks * theNThreads;
       
    //    NchunkG = (nPhotons/nDiv);
    //    if (NchunkG >= nthreads_total) {
    //       NchunkG -= ( (nPhotons/nDiv) % nthreads_total); // align on multiple of the grid size
    //    } else {
    //       fprintf(stderr,"Warning: too many division for photons (%d) resulting in low chunk size %d (< %d)\n",nDiv,NchunkG,nthreads_total);
    //    }

    //    NchunkE = (nElectrons/nDiv);
    //    if (NchunkE >= nthreads_total) {
    //       NchunkE -= ( (nElectrons/nDiv) % nthreads_total); // align on multiple of the grid size
    //    } else {
    //       fprintf(stderr,"Warning: too many division for electrons (%d) resulting in low chunk size %d (< %d)\n",nDiv,NchunkE,nthreads_total);
    //    }
       
    //    // cudaMalloc((void**)&photon_d0, NchunkG*sizeof(GXTrack));
    //    // cudaMalloc((void**)&photon_d1, NchunkG*sizeof(GXTrack));
       
    //    // cudaMalloc((void**)&electron_d0, NchunkE*sizeof(GXTrack));
    //    // cudaMalloc((void**)&electron_d1, NchunkE*sizeof(GXTrack));
    // }
    
    //stop time for cudaMalloc
    
    // cudaEventRecord (stop,0);
    // cudaEventSynchronize (stop); 
    // cudaEventElapsedTime (&elapsedTimeAlloc,start,stop);
    
    //start steam loop
    
    //prepare random engines on the device
    
    // curandState* devStates = 0;
    // cudaMalloc(&devStates, theNBlocks*theNThreads* sizeof(curandState));
    // curand_setup_gpu(devStates, time(NULL), theNBlocks, theNThreads);
    
    // curandState* devStates1 = 0;
    // cudaMalloc(&devStates1, theNBlocks*theNThreads* sizeof(curandState));
    // curand_setup_gpu(devStates1, time(NULL), theNBlocks, theNThreads);
    
    //start time for kernel
    // cudaEventRecord (start,0);
    
    //<<<---------------------------------------------------------->>>
    // if(kernelType ==1) {
    //   for(int i = 0, round = 0 ; i < nTracks ; i += 2*Nchunk, ++round) { 
        
    //     int count1 = min(Nchunk,nTracks-i);
    //     // fprintf(stderr,"round = %d count1 = %d count2 = %d\n",round,count1,-1);
    //     cudaMemcpyAsync(track_d0, track_h+i, count1*sizeof(GXTrack), 
    //                     cudaMemcpyHostToDevice, stream0);
        
    //     tracking_gpu(devStates,(GPVPhysicalVolume*)geom_d,bmap_d,
    //                  track_d0,
    //                  eBrem_table_d, eIoni_table_d, msc_table_d,
    //                  count1,
    //                  theNBlocks,theNThreads,stream0);
        
    //     cudaMemcpyAsync(track_h+i, track_d0, count1*sizeof(GXTrack), 
    //                     cudaMemcpyDeviceToHost, stream0);
        
    //     int count2 = min( Nchunk, nTracks - (i+count1) );
    //     if (count2) {
    //       cudaMemcpyAsync(track_d1, track_h+i+count1,count2*sizeof(GXTrack), 
    //                     cudaMemcpyHostToDevice, stream1);
        
    //       tracking_gpu(devStates,(GPVPhysicalVolume*)geom_d,bmap_d,
    //                    track_d1,
    //                    eBrem_table_d, eIoni_table_d, msc_table_d,
    //                    count2,
    //                    theNBlocks,theNThreads,stream1);
          
    //       cudaMemcpyAsync(track_h+i+count1, track_d1, count2*sizeof(GXTrack), 
    //                       cudaMemcpyDeviceToHost, stream1);
    //     }
    //     // fprintf(stderr,"round = %d count1 = %d count2 = %d\n",round,count1,count2);
    //     // random_testing_gpu( devStates, results,theNBlocks,theNThreads,stream0);
        
    //   }
    // }
    // else if(kernelType ==2) {
      
    //   //photon kernel
    //   for(int i = 0 ; i < nPhotons ; i += 2*NchunkG) { 
    //      int count1 = min(NchunkG,nPhotons-i) ;
    //      cudaMemcpyAsync(photon_d0, photon_h+i, count1*sizeof(GXTrack), 
    //                      cudaMemcpyHostToDevice, stream0);
         
    //      tracking_photon_gpu(devStates,(GPVPhysicalVolume*)geom_d,bmap_d,photon_d0,
    //                          eBrem_table_d, eIoni_table_d, msc_table_d,
    //                          count1, 
    //                          theNBlocks,theNThreads, stream0);
         
    //      cudaMemcpyAsync(photon_h+i, photon_d0, count1*sizeof(GXTrack), 
    //                      cudaMemcpyDeviceToHost, stream0);
         
    //      cudaMemcpyAsync(photon_d1, photon_h+i+NchunkG, (nPhotons/nDiv)*sizeof(GXTrack), 
    //                      cudaMemcpyHostToDevice, stream1);
         
    //      tracking_photon_gpu(devStates,(GPVPhysicalVolume*)geom_d,bmap_d,photon_d1,
    //                          eBrem_table_d, eIoni_table_d, msc_table_d,
    //                          (nPhotons/nDiv), theNBlocks,theNThreads, stream1);
         
    //      cudaMemcpyAsync(photon_h+i+NchunkG, photon_d1, (nPhotons/nDiv)*sizeof(GXTrack), 
    //                      cudaMemcpyDeviceToHost, stream1);
    //   }
      
    //   //electon kernel
    //   for(int i = 0 ; i < nElectrons ; i += 2*NchunkE) { 
    //      int count1 = min(NchunkE,nElectrons-i) ;

    //     cudaMemcpyAsync(electron_d0, electron_h+i, count1*sizeof(GXTrack), 
    //                     cudaMemcpyHostToDevice, stream0);
        
    //     tracking_electron_gpu(devStates1,(GPVPhysicalVolume*)geom_d,bmap_d,electron_d0,
    //                           eBrem_table_d, eIoni_table_d, msc_table_d,
    //                           count1, 
    //                           theNBlocks,theNThreads, stream0);
        
    //     cudaMemcpyAsync(electron_h+i, electron_d0, count1*sizeof(GXTrack), 
    //                     cudaMemcpyDeviceToHost, stream0);
        
    //     int count2 = min(NchunkE, nElectrons - (i+count1));
    //     cudaMemcpyAsync(electron_d1, electron_h+i+count1, count2*sizeof(GXTrack), 
    //                     cudaMemcpyHostToDevice, stream1);
        
    //     tracking_electron_gpu(devStates1,(GPVPhysicalVolume*)geom_d,bmap_d,electron_d1,
    //                           eBrem_table_d, eIoni_table_d, msc_table_d,
    //                           count2, 
    //                           theNBlocks,theNThreads, stream1);
        
    //     cudaMemcpyAsync(electron_h+i+count1, electron_d1, count2*sizeof(GXTrack), 
    //                     cudaMemcpyDeviceToHost, stream1);
    //   }
    // }
    
    //<<<---------------------------------------------------------->>>
    // cudaStreamSynchronize(stream0);
    // cudaStreamSynchronize(stream1);
    
    //stop time for kernel
    // cudaEventRecord (stop,0);
    // cudaEventSynchronize (stop);
    // cudaEventElapsedTime (&elapsedTimeGPU,start,stop);
    
    //start time for Memcpy (D2H)
    // cudaEventRecord (start,0);
    StopWatch timer;
    timer.start();

    //<<<---------------------------------------------------------->>>
    if(kernelType ==1) {
       tracking_cpu(geom->getWorldVolume(),*aNavigator,bmap,track_c,
                   eBrem_table,eIoni_table,msc_table,
                   nTracks);
    }
    else if(kernelType ==2 ) {
      tracking_photon_cpu(geom->getWorldVolume(),*aNavigator,bmap,
                          photon_c,
                          eBrem_table,eIoni_table,msc_table,
                          nPhotons);
      tracking_electron_cpu(geom->getWorldVolume(),*aNavigator,bmap,
                            electron_c,
                            eBrem_table,eIoni_table,msc_table,
                            nElectrons);
    }
    //<<<---------------------------------------------------------->>>

    timer.stop();
    // cudaEventRecord (stop,0);
    // cudaEventSynchronize (stop);
    
    float elapsedTime2 = timer.getTime();;
    // cudaEventElapsedTime (&elapsedTime2,start,stop);
    printf("Time Elapsed on CPU : %6.3f ms\n",elapsedTime2);
    
    //    printf("Ratio of Time Elapsed on CPU/GPU : %5.2f \n",
    //	   elapsedTime2/elapsedTimeGPU);
    // printf("Ratio of Time Elapsed on CPU/GPU : %5.2f \n",
    //        elapsedTime2/elapsedTimeGPU);
    
    //Compare results from GPU and CPU
    /*
     long int miscount = 0;
     for (size_t i = 0 ; i < nTracks ; i++) {
     if( track_h[i].x != track_c[i].x || track_h[i].y != track_c[i].y || 
     track_h[i].z != track_c[i].z ) {
     std::cout << "Compare results: check track position from gpu cpu " 
     << "(" << track_h[i].x << "," << track_h[i].y << "," 
     <<track_h[i].z << ")" << " (" << track_c[i].x << "," 
     << track_c[i].y <<","<<track_c[i].z << ")" << std::endl;
     miscount++;
     }
     }
     
     printf("Number of Cuda MisMatched Tracks = %ld\n",miscount);
     */
    //clean up: destory cuda event and free memory on device and host
    // cudaEventDestroy(start);
    // cudaEventDestroy(stop);
    
    if(kernelType == 1) {
      // cudaFree(track_d0);
      // cudaFree(track_d1);      
    }
    else if(kernelType ==2 || kernelType ==3) {
      // cudaFree(photon_d0);
      // cudaFree(photon_d1);
      
      // cudaFree(electron_d0);
      // cudaFree(electron_d1);
       
    }
    
    // cudaFree(devStates);
    // cudaFree(devStates1);
     delete [] steparr;
  }
  
  
  if(kernelType == 1) {
    // cudaFreeHost(track_h);
    // cudaFreeHost(track_c);
  }
  else if(kernelType ==2 || kernelType ==3) {
     // cudaFreeHost(photon_h);
     // cudaFreeHost(photon_c);
     // cudaFreeHost(electron_h);
     // cudaFreeHost(electron_c);
    // nothing
  }
  
  //end of event-loop
  // cudaStreamDestroy(stream0);
  // cudaStreamDestroy(stream1);
  
  //clean up: free memory on device and host
  // cudaFree(eBrem_table_d);
  // cudaFree(eIoni_table_d);
  // cudaFree(msc_table_d);
  
  // cudaFree(bmap_d);
  // free(bmap_h);
   
}
