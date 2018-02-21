/*
 mpic++ CMSApp_MPI.cc `root-config --cflags --glibs` -I/home/mba/GeantProject/GeantV/xsec/inc -I/home/mba/GeantProject/GeantV/vecprot_v2/inc/ -I/home/mba/GeantProject/GeantV/base/inc/ -I/home/mba/GeantProject/VecGeom/ -I/home/mba/GeantProject/VecGeom/VecCore/inc/ -I/home/mba/GeantProject/GeantV/vecprot_v2/src/ -I/home/mba/GeantProject/GeantV/examples/inc/ -I/home/mba/GeantV_install/hepmc3/include -L/home/mba/GeantProject/root-v06-03-02/lib/  -L/home/mba/GeantV_install/VecGeom/lib/ -L/home/mba/geant4-10.1/lib64/ -L/home/mba/GeantV_install/hepmc3/lib/ -L/home/mba/GeantV_install/GeantV/lib/ -lPhysics -lHist -lThread -lGeom -lVMC -lGeant_v -lXsec -lGeantExamples -lHepMC -lHepMCrootIO -o CMSApp_MPI
 
 NB: Event input file download here: https://www.dropbox.com/s/7ilyob22vpluy4c/pp14TEvminbias1000.root?dl=0
 */

#ifndef COPROCESSOR_REQUEST
#define COPROCESSOR_REQUEST false
#endif

#include <err.h>
#include <getopt.h>
#include <unistd.h>
#include "mpi.h"

#include "Rtypes.h"
#include "TGeoManager.h"

#include "GunGenerator.h"
#include "HepMCGenerator.h"
#include "TaskBroker.h"
#include "WorkloadManager.h"
#include "GeantPropagator.h"
#include "TTabPhysProcess.h"
#include "CMSApplication.h"
#include "HepMC/WriterRoot.h"
#include "HepMC/ReaderRoot.h"
#include "HepMC/GenEvent.h"
//#include "HepMC/HepMC.h"
#include "HepMCGenerator.h"

using namespace HepMC;

static int n_events = 10;
static int n_buffered = 5;
static int n_threads = 4;
static int n_track_max = 64;
static int n_learn_steps = 100000;
static int max_memory = 4000; /* MB */
static bool monitor = false, score = false, debug = false, coprocessor = false;

static struct option options[] = {{"events", required_argument, 0, 'e'},
                                  {"hepmc-event-file", required_argument, 0, 'E'},
                                  {"fstate", required_argument, 0, 'f'},
                                  {"geometry", required_argument, 0, 'g'},
                                  {"learn-steps", required_argument, 0, 'l'},
                                  {"max-tracks-per-basket", required_argument, 0, 'B'},
                                  {"monitor", no_argument, 0, 'm'},
                                  {"debug", no_argument, 0, 'd'},
                                  {"max-memory", required_argument, 0, 'M'},
                                  {"nbuffer", required_argument, 0, 'b'},
                                  {"score", no_argument, 0, 's'},
                                  {"threads", required_argument, 0, 't'},
                                  {"xsec", required_argument, 0, 'x'},
                                  {"coprocessor", required_argument, 0, 'r'},
                                  {0, 0, 0, 0}};

void help() {
  printf("\nUsage: cmsapp [OPTIONS] INPUT_FILE\n\n");

  for (int i = 0; options[i].name != NULL; i++) {
    printf("\t-%c  --%s\t%s\n", options[i].val, options[i].name, options[i].has_arg ? options[i].name : "");
  }
  printf("\n\n");
}


//______________________________________________________________________________
void readInputFile(int myRank, int nProc, int nEv, char *inputFile)
{
    //Es: nEv=32 total_events=100 nProc=2
    
    std::cout<<"readInputFile, my rank is "<<myRank<<"\n";
    //HepMC::ReaderRoot input_parser (inputFile);
    TString *fileNames[nProc];
    HepMC::WriterRoot * root_output[nProc];
    int events_parsed = 0;
    //int total_events = 0;
    /*
    while( !input_parser.failed() ) {
        HepMC::GenEvent evt(Units::GEV,Units::MM);
        input_parser.read_event(evt);
        if( input_parser.failed() ) break;
        total_events++;
    }
    input_parser.close();
    std::cout<<"readInputFile: total_events: "<<total_events<<"\n";*/
    
    HepMC::ReaderRoot text_input (inputFile);
    for (int i=0; i<nProc; i++)
    {
        fileNames[i]=new TString(Form("pp14TeVminbias1000_%d.root", i));
        root_output[i]=new HepMC::WriterRoot(fileNames[i]->Data());
    }
    
    //int nEv_perFile=total_events/nProc;
    int nEv_perFile=nEv/nProc; //nEv total divided by all the processes
    std::cout<<"We all need to process TOTAL nEvents "<<nEv<<" , and nEv per process "<<nEv_perFile<<"\n";
    //int resto=total_events%nProc;
    int index;
    //to debug
    int counter[nProc];
    for (int i=0; i<nProc; i++)
        counter[i]=0;
    
    while( !text_input.failed() ) {
        HepMC::GenEvent evt(Units::GEV,Units::MM);
        text_input.read_event(evt);
        if( text_input.failed() ) break;
        //if( events_parsed == 0 ) {
        //cout << "First event: " << endl;
        //Print::listing(evt);
        //}
        //++events_parsed;
        if(events_parsed<nEv)
        {
            index=events_parsed/nEv_perFile;
            if(index==nProc) index=nProc-1;
            //std::cout<<"events_parsed "<<events_parsed<<" index "<<index<<" nEv_perFile "<<nEv_perFile<<std::endl;
            counter[index]++;
            root_output[index]->write_event(evt);
            ++events_parsed;
        //if( events_parsed%1000 == 0 ) {
        //    cout << "Event: " << events_parsed << endl;
        //}
        } else break;
    }
    text_input.close();
    for(int i=0; i<nProc; i++)
    {
        root_output[i]->close();
        std::cout<<"counter["<<i<<"] "<<counter[i]<<std::endl;
    }
    std::cout << "Events parsed and written: " << events_parsed << "nEv: "<<nEv<<" nProc: "<<nProc<<" file splitted with nEv_perFile: "<<nEv_perFile<< std::endl;
    
}


//______________________________________________________________________________
void splitInputFile(int nProc, char *inputFile)
{
    //to count the total number of events
    //ReaderRoot input_parser ("pp14TeVminbias.root");
    
    std::cout<<"splitInputFile: Executing splitInputFile method.\n";
    HepMC::ReaderRoot input_parser (inputFile);
    TString *fileNames[nProc];
    HepMC::WriterRoot * root_output[nProc];
    int events_parsed = 0;
    int total_events = 0;
    while( !input_parser.failed() ) {
        HepMC::GenEvent evt(Units::GEV,Units::MM);
        input_parser.read_event(evt);
        if( input_parser.failed() ) break;
        total_events++;
    }
    input_parser.close();
    std::cout<<"splitInputFile: total_events: "<<total_events<<"\n";
    
    HepMC::ReaderRoot text_input (inputFile);
    for (int i=0; i<nProc; i++)
    {
        fileNames[i]=new TString(Form("pp14TeVminbias_%d.root", i));
        root_output[i]=new HepMC::WriterRoot(fileNames[i]->Data());
    }
    
    int nEv_perFile=total_events/nProc;
    //int resto=total_events%nProc;
    int index;
    //to debug
    //int counter[10];
    //for (int i=0; i<10; i++)
    //    counter[i]=0;
    
    while( !text_input.failed() ) {
        HepMC::GenEvent evt(Units::GEV,Units::MM);
        text_input.read_event(evt);
        if( text_input.failed() ) break;
        //if( events_parsed == 0 ) {
        //cout << "First event: " << endl;
        //Print::listing(evt);
        //}
        
        index=events_parsed/nEv_perFile;
        if(index==nProc) index=nProc-1;
        //std::cout<<"events_parsed "<<events_parsed<<" index "<<index<<" nEv_perFile "<<nEv_perFile<<std::endl;
        //counter[index]++;
        root_output[index]->write_event(evt);
        ++events_parsed;
        //if( events_parsed%1000 == 0 ) {
        //    cout << "Event: " << events_parsed << endl;
        //}
    }
    text_input.close();
    for(int i=0; i<nProc; i++)
    {
        root_output[i]->close();
        //std::cout<<"counter["<<i<<"] "<<counter[i]<<std::endl;
    }
    std::cout << "splitInputFile: Events parsed and written: " << events_parsed << std::endl;
    
}


int main(int argc, char *argv[]) {
  
    
    int nProc, myrank, namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    
    double starttime, endtime;
    MPI_Init(&argc, &argv);
    //starttime = MPI_Wtime();
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Get_processor_name(processor_name,&namelen);

    printf("From process %d out of %d, Hello World! I am running on processor %s\n", myrank, nProc, processor_name);
    system("pwd");
  

    std::string cms_geometry_filename("/home/mba/cms2015.root"); //absolute path of the file on the MPI server
    std::string xsec_filename("/home/mba/xsec_FTFP_BERT_G496p02_1mev.root"); //absolute path of the file on the MPI server
    std::string fstate_filename("/home/mba/fstate_FTFP_BERT_G496p02_1mev.root"); //absolute path of the file on the MPI server
    //std::string hepmc_event_filename(argv[1]);
    
    
    int nTOT_events=512;
    n_events=nTOT_events/nProc;
    //splitInputFile(nProc, argv[1]); //nProc input files, with n_TOT/nProc elements
    readInputFile(myrank, nProc, nTOT_events, argv[1]); //every process reads the input file, divides it in a number of files equal to the number of total processes
    std::string hepmc_event_filename(Form("pp14TeVminbias1000_%d.root", myrank)); //The i-th process elaborate its portion of the original file.
    
    //start measuring time
    starttime = MPI_Wtime();

    if (argc == 1) {
        help();
        exit(0);
    }

  while (true) {
    int c, optidx = 0;

    c = getopt_long(argc, argv, "E:e:f:g:l:B:mM:b:t:x:r", options, &optidx);

    if (c == -1)
      break;

    switch (c) {
    case 0:
      c = options[optidx].val;
    /* fall through */
    case 'e':
      n_events = (int)strtol(optarg, NULL, 10);

      if (n_events <= 0)
        errx(1, "number of events must be positive");
      break;

    case 'E':
      hepmc_event_filename = optarg;
      break;

    case 'f':
      fstate_filename = optarg;
      break;

    case 'g':
      cms_geometry_filename = optarg;
      break;

    case 'l':
      n_learn_steps = (int)strtol(optarg, NULL, 10);

      if (n_learn_steps <= 0)
        errx(1, "number of learning steps must be positive");
      break;

    case 'B':
      n_track_max = (int)strtol(optarg, NULL, 10);

      if (n_track_max < 1)
        errx(1, "max number of tracks per basket must be positive");
      break;

    case 'm':
      monitor = true;
      break;

    case 'M':
      max_memory = (int)strtol(optarg, NULL, 10);

      if (max_memory < 128)
        errx(1, "max memory is too low");
      break;

    case 'b':
      n_buffered = (int)strtol(optarg, NULL, 10);

      if (n_buffered < 1)
        errx(1, "number of buffered events must be positive");
      break;

    case 't':
      n_threads = (int)strtol(optarg, NULL, 10);

      if (n_threads < 1)
        errx(1, "number of threads must be positive");

      break;

    case 's':
      score = true;
      break;

    case 'x':
      xsec_filename = optarg;
      break;

    case 'r':
      coprocessor = optarg;
      break;

    default:
      errx(1, "unknown option %c", c);
    }
  }

  bool performance = true;
  TaskBroker *broker = nullptr;
  TGeoManager::Import(cms_geometry_filename.c_str());
  WorkloadManager *wmanager = WorkloadManager::Instance(n_threads);

  if (coprocessor) {
#ifdef GEANTCUDA_REPLACE
    CoprocessorBroker *gpuBroker = new CoprocessorBroker();
    gpuBroker->CudaSetup(32,128,1);
    broker = gpuBroker;
    nthreads += gpuBroker->GetNstream()+1;
#else
    std::cerr << "Error: Coprocessor processing requested but support was not enabled\n";
#endif
  }

  GeantPropagator *propagator = GeantPropagator::Instance(n_events, n_buffered);
  if (broker) propagator->SetTaskBroker(broker);
  wmanager->SetNminThreshold(5 * n_threads);
  propagator->fUseMonitoring = monitor;

  wmanager->SetMonitored(GeantPropagator::kMonQueue, monitor);
  wmanager->SetMonitored(GeantPropagator::kMonMemory, monitor);
  wmanager->SetMonitored(GeantPropagator::kMonBasketsPerVol, monitor);
  wmanager->SetMonitored(GeantPropagator::kMonVectors, monitor);
  wmanager->SetMonitored(GeantPropagator::kMonConcurrency, monitor);
  wmanager->SetMonitored(GeantPropagator::kMonTracksPerEvent, monitor);
  // Threshold for prioritizing events (tunable [0, 1], normally <0.1)
  // If set to 0 takes the default value of 0.01
  propagator->fPriorityThr = 0.1;

  // Initial vector size, this is no longer an important model parameter,
  // because is gets dynamically modified to accomodate the track flow
  propagator->fNperBasket = 16; // Initial vector size

  // This is now the most important parameter for memory considerations
  propagator->fMaxPerBasket = n_track_max;

  // Maximum user memory limit [MB]
  propagator->fMaxRes = max_memory;
  if (performance) propagator->fMaxRes = 0;
  propagator->fEmin = 0.001; // [10 MeV] energy cut
  propagator->fEmax = 0.01;  // 10 MeV

  propagator->LoadVecGeomGeometry();
  propagator->fProcess = new TTabPhysProcess("tab_phys", xsec_filename.c_str(), fstate_filename.c_str());

  if (hepmc_event_filename.empty()) {
    propagator->fPrimaryGenerator = new GunGenerator(propagator->fNaverage, 11, propagator->fEmax, -8, 0, 0, 1, 0, 0);
  } else {
    // propagator->fPrimaryGenerator->SetEtaRange(-2.,2.);
    // propagator->fPrimaryGenerator->SetMomRange(0.,0.5);
    // propagator->fPrimaryGenerator = new HepMCGenerator("pp14TeVminbias.hepmc3");
    propagator->fPrimaryGenerator = new HepMCGenerator(hepmc_event_filename);
  }
  propagator->fLearnSteps = n_learn_steps;
  if (performance) propagator->fLearnSteps = 0;

  CMSApplication *CMSApp = new CMSApplication();
  if (score) {
    CMSApp->SetScoreType(CMSApplication::kScore);
  } else {
    CMSApp->SetScoreType(CMSApplication::kNoScore);
  }
  propagator->fApplication = CMSApp;
  if (debug) {
    propagator->fUseDebug = true;
    propagator->fDebugTrk = 1;
  }
  propagator->fUseMonitoring = monitor;
  // Activate standard scoring   
  propagator->fUseStdScoring = true;
  if (performance) propagator->fUseStdScoring = false;
  propagator->PropagatorGeom(cms_geometry_filename.c_str(), n_threads, monitor);
    
    
    endtime = MPI_Wtime();
    std::cout<<"Execution Time: "<<endtime-starttime<<" sec."<<endl;
    MPI_Finalize();
  return 0;
}
