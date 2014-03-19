#include <iostream>
#include <fstream>
#include <iomanip>
#include "stdio.h"

#include <cuda.h>
#include <curand.h>
#include <cuda_runtime.h>

#include "GXTrack.h"
#include "GPEMPhysicsUtils.h"
#include "electronTest_kernel.h"
#include "random_kernel.h"
#include "dma_kernel.h"

#include "GPHistoManager.hh"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

#include "GPPhysics2DVector.h"

using namespace std;

int main (int argc, char* argv[])
{
  //argument: kernel type 
  int emModel = 0; //0 electron 1 brem, 2 ioni, 3 msc
                   //21 brem, 22 ioni to stack secondaries with a dynamic memory
  int runType = 0; //0 normal, 1 test
  int isStack = 0; //0 stack secondaries, 1 no secondary stack
  int trkReadMode = 0; //0 randomly generate track parameters,  1 read track parameters from file

  if(argc >= 2) emModel = atoi(argv[1]);
  if(argc >= 3) runType = atoi(argv[2]);
  if(argc >= 4) isStack = atoi(argv[3]);
  if(argc >= 5) trkReadMode = atoi(argv[4]);

  if(emModel < 0 || emModel > 40) {
    std::cout << "Usage: electronTest model=[0-3,21-22] [0-1] [0-1]: ... " 
	      << std::endl;
    return 0;
  }

  int nDevice;
  bool cudaEnabled = false;

  cudaGetDeviceCount(&nDevice);

  if(nDevice > 0) {
    cudaDeviceReset();
    cudaEnabled = true;
  }
  else {
    std::cout << "Warning: No Cuda Capable Device ... " << std::endl;
  }

  // set the default number of threads and thread blocks
  int theNBlocks  =  32;
  int theNThreads = 128;

  if(runType == 1 ) {
    theNBlocks  =  2;
    theNThreads =  4;
  }

  //@@@pre-allocated memory for Secondaries
  //maximum size of dynamic memory allocation on the device
  //set 256 megabytes on the heap (global memory)
  cudaThreadSetLimit(cudaLimitMallocHeapSize, 256*1024*1024);

  // prepare physics table
  bool useSpline = true;

  GPPhysicsTable eBrem_table;
  readTable(&eBrem_table,"data/Lambda.eBrem.e-.asc");
  int nv = eBrem_table.nPhysicsVector;
  for(int j=0; j < nv; j++){
    eBrem_table.physicsVectors[j].SetSpline(useSpline);
  }

  //Ionisation

  GPPhysicsTable eIoni_table;
  readTable(&eIoni_table,"data/Lambda.eIoni.e-.asc");
  nv = eIoni_table.nPhysicsVector;
  for(int j=0; j < nv; j++){
    eIoni_table.physicsVectors[j].SetSpline(useSpline);
  }

  GPPhysicsTable eIoni_range;
  readTable(&eIoni_range,"data/Range.eIoni.e-.asc");
  nv = eIoni_range.nPhysicsVector;
  for(int j=0; j < nv; j++){
    eIoni_range.physicsVectors[j].SetSpline(useSpline);
  }

  GPPhysicsTable eIoni_dedx;
  readTable(&eIoni_dedx,"data/DEDX.eIoni.e-.asc");
  nv = eIoni_dedx.nPhysicsVector;
  for(int j=0; j < nv; j++){
    eIoni_dedx.physicsVectors[j].SetSpline(useSpline);
  }

  GPPhysicsTable eIoni_invr;
  readTable(&eIoni_invr,"data/InverseRange.eIoni.e-.asc");
  nv = eIoni_invr.nPhysicsVector;
  for(int j=0; j < nv; j++){
    eIoni_invr.physicsVectors[j].SetSpline(useSpline);
  }

  //msc

  GPPhysicsTable msc_table;
  readTable(&msc_table,"data/Lambda.msc.e-.asc");
  nv = msc_table.nPhysicsVector;
  for(int j=0; j < nv; j++){
    msc_table.physicsVectors[j].SetSpline(useSpline);
  }

  printf("Copying Physics data tables from host to device\n");

  GPPhysicsTable* eBrem_table_d;

  GPPhysicsTable* eIoni_table_d;
  GPPhysicsTable* eIoni_range_d;
  GPPhysicsTable* eIoni_dedx_d;
  GPPhysicsTable* eIoni_invr_d;

  GPPhysicsTable* msc_table_d;

  if(cudaEnabled) {
    cudaMalloc((void**)&eBrem_table_d, sizeof(GPPhysicsTable));
    cudaMemcpy(eBrem_table_d, &eBrem_table, sizeof(GPPhysicsTable),
	       cudaMemcpyHostToDevice);

    cudaMalloc((void**)&eIoni_table_d, sizeof(GPPhysicsTable));
    cudaMemcpy(eIoni_table_d, &eIoni_table, sizeof(GPPhysicsTable),
	       cudaMemcpyHostToDevice);

    cudaMalloc((void**)&eIoni_range_d, sizeof(GPPhysicsTable));
    cudaMemcpy(eIoni_range_d, &eIoni_range, sizeof(GPPhysicsTable),
	       cudaMemcpyHostToDevice);

    cudaMalloc((void**)&eIoni_dedx_d, sizeof(GPPhysicsTable));
    cudaMemcpy(eIoni_dedx_d, &eIoni_dedx, sizeof(GPPhysicsTable),
	       cudaMemcpyHostToDevice);

    cudaMalloc((void**)&eIoni_invr_d, sizeof(GPPhysicsTable));
    cudaMemcpy(eIoni_invr_d, &eIoni_invr, sizeof(GPPhysicsTable),
	       cudaMemcpyHostToDevice);

    cudaMalloc((void**)&msc_table_d, sizeof(GPPhysicsTable));
    cudaMemcpy(msc_table_d, &msc_table, sizeof(GPPhysicsTable),
	       cudaMemcpyHostToDevice);
  }
  
  //G4SeltzerBergerModel data
  G4int maxZ = 92;
  GPPhysics2DVector* sbData =
    (GPPhysics2DVector*) malloc(maxZ*sizeof(GPPhysics2DVector));

  char sbDataFile[256];
  for(G4int iZ = 0 ; iZ < maxZ ; iZ++) {  
    sprintf(sbDataFile,"brem_SB/br%d",iZ+1);
    std::ifstream fin(sbDataFile);
    G4bool check = RetrieveSeltzerBergerData(fin, &sbData[iZ]);
    if(!check) {
      printf("Failed To open SeltzerBerger Data file for Z= %d\n",iZ+1);
    }
  }

  GPPhysics2DVector* sbData_d;
  if(cudaEnabled) {
    cudaMalloc((void**)&sbData_d,maxZ*sizeof(GPPhysics2DVector));
    cudaMemcpy(sbData_d, sbData, maxZ*sizeof(GPPhysics2DVector),
	       cudaMemcpyHostToDevice);
  }

  //=======  Booking histograms =================
  GPHistoManager& hmgr = GPHistoManager::getInstance();

  TH1F& hdp  = hmgr.book1F("hdp", 100,0,1000);
  TH1F& hdpx = hmgr.book1F("hdpx",100,-1200,1200);
  TH1F& hdpy = hmgr.book1F("hdpy",100,-1200,1200);
  TH1F& hdpz = hmgr.book1F("hdpz",100,-1000,1000);
  TH1F& hdE  = hmgr.book1F("hdE", 100,0,1200);
  TH1F& hdq  = hmgr.book1F("hdq", 5, -2.5, 2.5);
  TH1F& hds  = hmgr.book1F("hds",100, 0, 0.5);
  TH1F& hdx  = hmgr.book1F("hdx",100,-1600,1600);
  TH1F& hdy  = hmgr.book1F("hdy",100,-1600,1600);
  TH1F& hdz  = hmgr.book1F("hdz",100,-3200,3200);
  TH1F& hdnsec = hmgr.book1F("hdnsec",100,0,30000);
  TH2F& hdxy = hmgr.book2F("hdxy",100,-1600,1600,100,-1600,1600);

  TH1F& hcp  = hmgr.book1F("hcp",100,0,1000);
  TH1F& hcpx = hmgr.book1F("hcpx",100,-1200,1200);
  TH1F& hcpy = hmgr.book1F("hcpy",100,-1200,1200);
  TH1F& hcpz = hmgr.book1F("hcpz",100,-1000,1000);
  TH1F& hcE  = hmgr.book1F("hcE",100,0,1200);
  TH1F& hcq  = hmgr.book1F("hcq", 5, -2.5, 2.5);
  TH1F& hcs  = hmgr.book1F("hcs",100, 0, 0.5);
  TH1F& hcx  = hmgr.book1F("hcx",100,-1600,1600);
  TH1F& hcy  = hmgr.book1F("hcy",100,-1600,1600);
  TH1F& hcz  = hmgr.book1F("hcz",100,-3200,3200);
  TH1F& hcnsec = hmgr.book1F("hcnsec",100,0,30000);
  TH2F& hcxy = hmgr.book2F("hcxy",100,-1600,1600,100,-1600,1600);

  TH1F& hsdp  = hmgr.book1F("hsdp",100,0,1000);
  TH1F& hsdpx = hmgr.book1F("hsdpx",100,-1200,1200);
  TH1F& hsdpy = hmgr.book1F("hsdpy",100,-1200,1200);
  TH1F& hsdpz = hmgr.book1F("hsdpz",100,-1000,1000);
  TH1F& hsdE  = hmgr.book1F("hsdE",100,0,1200);
  TH1F& hsdq  = hmgr.book1F("hsdq", 5, -2.5, 2.5);
  TH1F& hsds  = hmgr.book1F("hsds",100, 0, 0.5);
  TH1F& hsdx  = hmgr.book1F("hsdx",100,-1600,1600);
  TH1F& hsdy  = hmgr.book1F("hsdy",100,-1600,1600);
  TH1F& hsdz  = hmgr.book1F("hsdz",100,-3200,3200);
  TH2F& hsdxy = hmgr.book2F("hsdxy",100,-1600,1600,100,-1600,1600);

  TH1F& hscp  = hmgr.book1F("hscp",100,0,1000);
  TH1F& hscpx = hmgr.book1F("hscpx",100,-1200,1200);
  TH1F& hscpy = hmgr.book1F("hscpy",100,-1200,1200);
  TH1F& hscpz = hmgr.book1F("hscpz",100,-1000,1000);
  TH1F& hscE  = hmgr.book1F("hscE",100,0,1200);
  TH1F& hscq  = hmgr.book1F("hscq", 5, -2.5, 2.5);
  TH1F& hscs  = hmgr.book1F("hscs",100, 0, 0.5);
  TH1F& hscx  = hmgr.book1F("hscx",100,-1600,1600);
  TH1F& hscy  = hmgr.book1F("hscy",100,-1600,1600);
  TH1F& hscz  = hmgr.book1F("hscz",100,-3200,3200);
  TH2F& hscxy = hmgr.book2F("hscxy",100,-1600,1600,100,-1600,1600);

  // book histograms for CPU (host kernel)
  hmgr.book1F("hioniStepLength",100,0,10);
  hmgr.book1F("hstepDiff",100,0,1);

  if(emModel==1) hmgr.bookBremHistos();
  if(emModel==2) hmgr.bookIoniHistos();
  if(emModel==3) hmgr.bookMscHistos();

  //=======  end of histogram booking ============

  // get the list of tracks from the stack
  printf("Preparing input tracks at host\n");

  //event-loop
  int nevent = 10;
  ifstream* trkfile = 0;
  if(trkReadMode) {
    trkfile = new ifstream("trkfile");
    string line;
    std::getline(*trkfile, line);  // skip first line
  }

  for(size_t i = 0; i<10; ++i) {
    cout<<" Checking Rndm(): "<< Rndm() << endl;
  }

  for (int i = 0 ; i < nevent ; i++) {

    int nTracks = 20000;
    if(runType == 1 ) nTracks = 2;
#ifdef GPUDEBUG
    nTracks = 2;
#endif

    printf("Event %d> Transportation in kernel for %d tracks\n",i,nTracks);

    // populate tracks with track position, momentum, energy, steplength
    GXTrack *track_h;
    cudaHostAlloc((void**) &track_h, nTracks*sizeof(GXTrack),
                  cudaHostAllocDefault);

    //repeat on the CPU side, but with malloc
    GXTrack *track_c = (GXTrack *) malloc (nTracks*sizeof(GXTrack));

    //generate tracks with random position and momentum

    G4double photonFraction = 0.0;
    const G4double kinEmin = 20.;
    const G4double kinEmax = 1000.;

    G4double rho, z, p, mass, kinE;
    G4double theta, phi;

    for(size_t i = 0 ; i < nTracks ; i++){

      if(trkfile) {
	// read tracks from file
	(*trkfile) >> track_h[i].q >> track_h[i].x >> track_h[i].y >> track_h[i].z >> track_h[i].s
		   >> track_h[i].px >> track_h[i].py >> track_h[i].pz >> track_h[i].E;
	if(i<3) cout <<"Track read: "<< track_h[i].q <<"   "<< track_h[i].x <<"   "<< track_h[i].y <<"   "<< track_h[i].z <<"   "<< track_h[i].s
		     <<"   "<< track_h[i].px <<"   "<< track_h[i].py <<"   "<< track_h[i].pz <<"   "<< track_h[i].E << endl;
      }
      else {
	// random track parameters
	rho = ecalRmim + (ecalRmax-ecalRmim)*Rndm();
	z = ecalZmax*(2*Rndm()-1.0);
	phi = twopi*Rndm();
	theta = std::atan(rho/z);
	if(theta<0) theta += pi;

	track_h[i].status = 0;
	track_h[i].q = ( Rndm() < photonFraction ) ? 0.0 : -1.0;
	track_h[i].x = rho*std::cos(phi);
	track_h[i].y = rho*std::sin(phi); 
	track_h[i].z = z; 
	track_h[i].s = 10.0*Rndm(); 

	kinE = kinEmin + (kinEmax - kinEmin)*Rndm();
	mass = electron_mass_c2*track_h[i].q*track_h[i].q;
	p = sqrt(kinE*(kinE+2*mass));
	track_h[i].px = p*std::sin(theta)*std::cos(phi);
	track_h[i].py = p*std::sin(theta)*std::sin(phi);
	track_h[i].pz = p*std::cos(theta);

	//track_h[i].E  = p*p/(sqrt(p*p + mass*mass) + mass);
	//	track_h[i].E  = sqrt(p*p + mass*mass) - mass;  // this is simpler!
	track_h[i].E = kinE;
      }

      CopyTrack(&track_h[i],&track_c[i]);
    }

    float elapsedTimeAlloc = 0.0;
    float elapsedTimeUp = 0.0;
    float elapsedTimeKernel = 0.0;
    float elapsedTimeDown = 0.0;
    float elapsedTimeGPU = 0.0;
    
    printf("Copying track data from host to device\n");
    
    //start time memory allocation on GPU
    cudaEvent_t start, stop;
    cudaEventCreate (&start);
    cudaEventCreate (&stop);

    //start time for allocation
    cudaEventRecord (start,0);
    
    GXTrack *track_d;
    cudaMalloc((void**)&track_d, nTracks*sizeof(GXTrack));
    
    //stop time for cudaMalloc
    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);
    cudaEventElapsedTime (&elapsedTimeAlloc,start,stop);
    
    //start time for Memcpy (H2D)
    cudaEventRecord (start,0);

    cudaMemcpy(track_d, track_h, nTracks*sizeof(GXTrack), 
	       cudaMemcpyHostToDevice);
    
    //stop time for Memcpy
    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);
    cudaEventElapsedTime (&elapsedTimeUp,start,stop);

    //-- measure time to setup memory for array of secondaries + random engines
    float elapsedTimeSecondaries = 0.0;
    cudaEventRecord (start,0);

    //prepare random engines on the device
    curandState* devStates = 0;
    cudaMalloc(&devStates, theNBlocks*theNThreads* sizeof(curandState));
    curand_setup_gpu(devStates, time(NULL), theNBlocks, theNThreads);

    //@@@pre-allocated memory for secondaries
    G4int stackSize = 0;
    GXTrack *secTracks_d;

    cudaMalloc((void**)&secTracks_d, 
	       maxSecondaryPerStep*nTracks*sizeof(GXTrack));

    //atomic counter for the total number of secondaries
    G4int *stackSize_d;
    cudaMalloc((void**)&stackSize_d, sizeof(G4int));
    cudaMemset(stackSize_d,0,sizeof(G4int));

    //atomic counter for the last array position of secondaries
    G4int *offset_d;
    cudaMalloc((void**)&offset_d, sizeof(G4int));
    cudaMemset(offset_d,0,sizeof(G4int));

    //fixed size of memory for secondaries on host
    GXTrack *secTracks_d2h
      = (GXTrack*) malloc(maxSecondaryPerStep*nTracks*sizeof(GXTrack));

    GXTrack *secTracks_h
      = (GXTrack*) malloc(maxSecondaryPerStep*nTracks*sizeof(GXTrack));

    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);
    cudaEventElapsedTime (&elapsedTimeSecondaries,start,stop);



    //start time for kernel
    cudaEventRecord (start,0);
    //<<<---------------------------------------------------------->>>
    if(cudaEnabled) {
      if(emModel == 0 ) {

        electron_gpu(devStates, track_d, nTracks,
		     eBrem_table_d, eIoni_table_d,
		     eIoni_range_d, eIoni_dedx_d, eIoni_invr_d, 
		     msc_table_d, sbData_d,
		     secTracks_d,stackSize_d,offset_d,isStack,runType,
		     theNBlocks,theNThreads);
        if(isStack==0) {
          cudaThreadSynchronize();
          cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),
                     cudaMemcpyDeviceToHost);
	  printf("stackSize=%d\n", stackSize);
        }         
      }
      else if(emModel == 1 ) {

	brem_gpu(devStates, track_d, nTracks,
		 eBrem_table_d, eIoni_table_d, msc_table_d, sbData_d,
		 secTracks_d,stackSize_d,offset_d,isStack,runType,
		 theNBlocks,theNThreads);

        if(isStack==0) {
          cudaThreadSynchronize();
          cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),
                     cudaMemcpyDeviceToHost);
        }
      }
      else if(emModel == 2 ) {

	ioni_gpu(devStates, track_d, nTracks,
		 eBrem_table_d, eIoni_table_d, 
		 eIoni_range_d, eIoni_dedx_d, eIoni_invr_d, 
		 msc_table_d, sbData_d,
		 secTracks_d,stackSize_d,offset_d,isStack,runType,
		 theNBlocks,theNThreads);
        if(isStack==0) {
          cudaThreadSynchronize();
          cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),
                     cudaMemcpyDeviceToHost);
        }       
      }
      else if(emModel == 3 ) {

	msc_gpu(devStates, track_d, nTracks,
		eBrem_table_d, eIoni_table_d, msc_table_d, sbData_d,
		secTracks_d,stackSize_d,offset_d,isStack,runType,
		theNBlocks,theNThreads);
        if(isStack==0) {
          cudaThreadSynchronize();
          cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),
                     cudaMemcpyDeviceToHost);
        }         
      }
      //<<<---------------------------------------------------------->>>
      else if(emModel == 11 ) {

	//memory allocation for a containter holding secondaries per block 
	GXSecContainer *secContainer_d;
	cudaMalloc((void**)&secContainer_d, 
		   theNBlocks*sizeof(GXSecContainer));

	brem_gpu_dma(devStates, track_d, nTracks,
		     eBrem_table_d, eIoni_table_d, msc_table_d, sbData_d,
		     secContainer_d,stackSize_d,offset_d,
		     theNBlocks,theNThreads);

	cudaThreadSynchronize();

	//get the number of secondaries created on GPU
	cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),
		   cudaMemcpyDeviceToHost);

	cudaMalloc((void**)&secTracks_d, stackSize*sizeof(GXTrack));

	//reset the size of stackSize
	stackSize = 0;
	cudaMemset(stackSize_d,0,sizeof(G4int));

	//reallocate secondaries into one array
	//<---------------------------------------------------------------------
	dma_realloc_perblock_gpu(nTracks,secContainer_d,secTracks_d,stackSize_d,
				 offset_d,theNBlocks,1);
	//<---------------------------------------------------------------------
	cudaThreadSynchronize();

	//get the number of reallocated secondaries on GPU
	cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),cudaMemcpyDeviceToHost);
	dma_free_perblock_gpu(nTracks,secContainer_d,theNBlocks,1);
	cudaFree(secContainer_d);

      }
      else if(emModel == 12 ) {

	//memory allocation for a containter holding secondaries per block 
	GXSecContainer *secContainer_d;
	cudaMalloc((void**)&secContainer_d, 
		   theNBlocks*sizeof(GXSecContainer));

	ioni_gpu_dma(devStates, track_d, nTracks,
		     eBrem_table_d, eIoni_table_d, 
		     eIoni_range_d, eIoni_dedx_d, eIoni_invr_d, 
		     msc_table_d, sbData_d,
		     secContainer_d,stackSize_d,offset_d,
		     theNBlocks,theNThreads);

	cudaThreadSynchronize();

	//get the number of secondaries created on GPU
	cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),
		   cudaMemcpyDeviceToHost);

	cudaMalloc((void**)&secTracks_d, stackSize*sizeof(GXTrack));

	//reset the size of stackSize
	stackSize = 0;
	cudaMemset(stackSize_d,0,sizeof(G4int));

	//reallocate secondaries into one array
	//<---------------------------------------------------------------------
	dma_realloc_perblock_gpu(nTracks,secContainer_d,secTracks_d,stackSize_d,
				 offset_d,theNBlocks,1);
	//<---------------------------------------------------------------------
	cudaThreadSynchronize();

	//get the number of reallocated secondaries on GPU
	cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),cudaMemcpyDeviceToHost);
	dma_free_perblock_gpu(nTracks,secContainer_d,theNBlocks,1);
	cudaFree(secContainer_d);
      }
    }

    cudaThreadSynchronize();
  
    //stop time for kernel
    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);
    cudaEventElapsedTime (&elapsedTimeKernel,start,stop);

    //<<<---------------------------------------------------------->>>

    // fill histograms
    G4double px,py,pz;
    for(size_t i = 0 ; i < nTracks ; i++){
      px = track_h[i].px;
      py = track_h[i].py;
      pz = track_h[i].pz;
      hdp.Fill( sqrt(px*px+py*py+pz*pz) ); 
      hdpx.Fill( px );
      hdpy.Fill( py );
      hdpz.Fill( pz );
      hdE.Fill( track_h[i].E );
      hdq.Fill( track_h[i].q );
      hds.Fill( track_h[i].s );
      hdx.Fill( track_h[i].x );
      hdy.Fill( track_h[i].y );
      hdz.Fill( track_h[i].z );
      hdxy.Fill(track_h[i].x,track_h[i].y);
      
      px = track_c[i].px;
      py = track_c[i].py;
      pz = track_c[i].pz;
      hcp.Fill( sqrt(px*px+py*py+pz*pz) ); 
      hcpx.Fill( px );
      hcpy.Fill( py );
      hcpz.Fill( pz );
      hcE.Fill( track_c[i].E );
      hcq.Fill( track_c[i].q );
      hcs.Fill( track_c[i].s );
      hcx.Fill( track_c[i].x );
      hcy.Fill( track_c[i].y );
      hcz.Fill( track_c[i].z );
      hcxy.Fill( track_c[i].x, track_c[i].y );
    }

    //<<<---------------------------------------------------------->>>
    
    //start time for Memcpy (D2H)
    cudaEventRecord (start,0);
    
    if(cudaEnabled) {
    	cudaMemcpy(track_h, track_d,nTracks*sizeof(GXTrack), 
		   cudaMemcpyDeviceToHost);

	// GL: experimenting to obtain the secondary tracks
	printf("stackSize=%d\n", stackSize);
    	cudaMemcpy(secTracks_d2h, secTracks_d,stackSize*sizeof(GXTrack), 
		   cudaMemcpyDeviceToHost);

	// for(int i=0; i<stackSize+3; ++i) {
	//   printf("i=%i, E=%10.6e, pz=%10.6e\n",i,secTracks_h[i].E, secTracks_h[i].pz);
	// }
    }

    //stop time for Memcpy
    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);
    cudaEventElapsedTime (&elapsedTimeDown,start,stop);

    elapsedTimeGPU = elapsedTimeUp + elapsedTimeKernel + elapsedTimeDown + elapsedTimeAlloc;
    
    printf("Time Elapsed on GPU : %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f ms\n",
	   elapsedTimeAlloc, elapsedTimeUp,elapsedTimeKernel,
	   elapsedTimeDown,elapsedTimeGPU,elapsedTimeSecondaries);

    //<<<---------------------------------------------------------->>>

    printf("Executing the host code in CPU\n");

    //reset stackSize for CPU
    G4int stackSize_c = 0;
   
    cudaEventRecord (start,0);
    
    if(emModel == 0 ) {
      electron_cpu(track_c,nTracks,
		   &eBrem_table,&eIoni_table,
		   &eIoni_range, &eIoni_dedx, &eIoni_invr,&msc_table,
		   sbData,secTracks_h,&stackSize_c,
		   isStack,runType);
    }
    else if(emModel == 1 ) {
      brem_cpu(track_c,nTracks,
	       &eBrem_table,&eIoni_table,&msc_table,sbData,
	       secTracks_h,&stackSize_c,
	       isStack,runType);
    }
    else if(emModel == 2 ) {
      ioni_cpu(track_c,nTracks,
	       &eBrem_table,&eIoni_table,
	       &eIoni_range, &eIoni_dedx, &eIoni_invr, &msc_table,
	       sbData,secTracks_h,&stackSize_c,
	       isStack,runType);
    }
    else if(emModel == 3 ) {
      msc_cpu(track_c,nTracks,
	      &eBrem_table,&eIoni_table,&msc_table,sbData,
	      secTracks_h,&stackSize_c,isStack,runType);
    }
    //<<<---------------------------------------------------------->>>
    else if(emModel == 11 ) {
      brem_cpu_dma(track_c,nTracks,
	       &eBrem_table,&eIoni_table,&msc_table,sbData);
    }
    else if(emModel == 12 ) {
      ioni_cpu_dma(track_c,nTracks,
	       &eBrem_table,&eIoni_table,
	       &eIoni_range, &eIoni_dedx, &eIoni_invr, 
	       &msc_table,sbData);
    }
    //<<<---------------------------------------------------------->>>

    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);

    float elapsedTime2 = 0.0;
    cudaEventElapsedTime (&elapsedTime2,start,stop);
    printf("Time Elapsed on CPU : %6.3f ms\n",elapsedTime2);

    printf("Ratio of Time Elapsed on CPU/GPU : %5.2f %5.2f \n",
	   elapsedTime2/elapsedTimeGPU,elapsedTime2/elapsedTimeKernel);


    // fill secTracks histograms
    hdnsec.Fill(stackSize);
    hcnsec.Fill(stackSize_c);
    for(size_t i = 0 ; i < stackSize ; i++){
      px = secTracks_d2h[i].px;
      py = secTracks_d2h[i].py;
      pz = secTracks_d2h[i].pz;
      hsdp.Fill( sqrt(px*px+py*py+pz*pz) ); 
      hsdpx.Fill( px );
      hsdpy.Fill( py );
      hsdpz.Fill( pz );
      hsdE.Fill( secTracks_d2h[i].E );
      hsdq.Fill( secTracks_d2h[i].q );
      hsds.Fill( secTracks_d2h[i].s );
      hsdx.Fill( secTracks_d2h[i].x );
      hsdy.Fill( secTracks_d2h[i].y );
      hsdz.Fill( secTracks_d2h[i].z );
      hsdxy.Fill( secTracks_d2h[i].x,secTracks_d2h[i].y);

      px = secTracks_h[i].px;
      py = secTracks_h[i].py;
      pz = secTracks_h[i].pz;
      hscp.Fill( sqrt(px*px+py*py+pz*pz) ); 
      hscpx.Fill( px );
      hscpy.Fill( py );
      hscpz.Fill( pz );
      hscE.Fill( secTracks_h[i].E );
      hscq.Fill( secTracks_h[i].q );
      hscs.Fill( secTracks_h[i].s );
      hscx.Fill( secTracks_h[i].x );
      hscy.Fill( secTracks_h[i].y );
      hscz.Fill( secTracks_h[i].z );
      hscxy.Fill( secTracks_h[i].x,secTracks_h[i].y);
    }

    //clean up: destroy cuda event and free memory on device and host
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    if(cudaEnabled) {
      cudaFree(track_d);
      cudaFree(devStates);
    }
    cudaFreeHost(track_h);
    free(track_c);
    free(secTracks_h);
    free(secTracks_d2h);

    cudaFree(stackSize_d);
    cudaFree(offset_d);
    cudaFree(secTracks_d);
  }

  //end of event-loop
  if(cudaEnabled) {
    cudaFree(eBrem_table_d);
    cudaFree(eIoni_table_d);
    cudaFree(eIoni_range_d);
    cudaFree(eIoni_dedx_d);
    cudaFree(eIoni_invr_d);
    cudaFree(msc_table_d);
  }

  free(sbData);
  cudaFree(sbData_d);

  // save histos
  //printf("Saving histograms...\n");
  //  hmgr.saveHistos();
  printf("Destroying HistoManager...\n");
  hmgr.destroy();
}
