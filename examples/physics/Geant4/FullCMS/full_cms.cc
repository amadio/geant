
//--------------------------------------------------------
// Geant4 full CMS application: 22 November 2017 (README)
//--------------------------------------------------------

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "G4UIsession.hh"
#include "G4UIterminal.hh"

#include "G4PhysListFactory.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VModularPhysicsList.hh"

#include "Randomize.hh"
#include "MyDetectorConstruction.hh"
#include "MyGVPhysicsList.hh"

#include "MyActionInitialization.hh"

#include <getopt.h>
#include <err.h>
#include <iostream>
#include <iomanip>

static bool parIsPerformance        = false;
static std::string parMacroFileName = "";
static std::string parPhysListName  = "GV";

void GetInputArguments(int argc, char **argv);
void Help();

int main(int argc, char **argv)
{
  //
  // get input arguments
  GetInputArguments(argc, argv);
  G4cout << " ========== Running full_cms ========================= " << G4endl
         << "   Physics list name   =  " << parPhysListName << G4endl
         << "   Geant4 macro        =  " << parMacroFileName << G4endl
         << "   Performance mode    =  " << parIsPerformance << G4endl
         << " ===================================================== " << G4endl;
  //
  // choose the Random engine: set to MixMax explicitely (default form 10.4)
  G4Random::setTheEngine(new CLHEP::MixMaxRng);
  // set seed and print info
  G4Random::setTheSeed(12345678);
  G4cout << G4endl << " ===================================================== " << G4endl
         << " Random engine = " << G4Random::getTheEngine()->name() << G4endl
         << " Initial seed  = " << G4Random::getTheSeed() << G4endl
         << " ===================================================== " << G4endl << G4endl;
//
// Construct the default run manager
#ifdef G4MULTITHREADED
  G4MTRunManager *runManager = new G4MTRunManager;
  // Number of threads can be defined via macro command
  G4int nThreads = 4;
  runManager->SetNumberOfThreads(nThreads);
#else
  G4RunManager *runManager = new G4RunManager;
#endif
  //
  // set mandatory initialization classes
  //
  // 1. Detector construction
  MyDetectorConstruction *detector = new MyDetectorConstruction;
  runManager->SetUserInitialization(detector);
  //
  // 2. Physics list
  G4PhysListFactory factory;
  if (factory.IsReferencePhysList(parPhysListName)) {
    G4VModularPhysicsList *physList = factory.GetReferencePhysList(parPhysListName);
    runManager->SetUserInitialization(physList);
  } else if (parPhysListName == G4String("GV")) {
    G4VUserPhysicsList *physList = new MyGVPhysicsList();
    runManager->SetUserInitialization(physList);
  } else {
    G4cerr << "ERROR: Physics List " << parPhysListName << " UNKNOWN!" << G4endl;
    return -1;
  }
  //
  // 3. User action
  runManager->SetUserInitialization(new MyActionInitialization(parIsPerformance));
  //
  // 4. Run the simulation in batch mode
  G4UImanager *UI  = G4UImanager::GetUIpointer();
  G4String command = "/control/execute ";
  UI->ApplyCommand(command + parMacroFileName);
  //
  // Print out the final random number
  G4cout << G4endl << " ================================================================= " << G4endl
         << " Final random number = " << CLHEP::HepRandom::getTheEngine()->flat() << G4endl
         << " ================================================================= " << G4endl << G4endl;
  //
  // Delete the RunManager
  delete runManager;
  return 0;
}

static struct option options[] = {
    {"standard Geant4 macro file (REQUIRED)", required_argument, 0, 'm'},
    {"physics list name (default GV)", required_argument, 0, 'f'},
    {"flag to run the application in performance mode (default FALSE)", no_argument, 0, 'p'},
    {0, 0, 0, 0}};

void Help()
{
  std::cout << "\n " << std::setw(100) << std::setfill('=') << "" << std::setfill(' ') << std::endl;
  std::cout << "  Full CMS Geant4 application.    \n"
            << std::endl
            << "  **** Parameters: \n"
            << "      -m :   REQUIRED : the standard Geamt4 macro file \n"
            << "      -f :   physics list name (default: GV) \n"
            << "      -p :   flag  ==> run the application in performance mode i.e. no user actions \n"
            << "         :   -     ==> run the application in NON performance mode i.e. with user actions (default) \n"
            << std::endl;
  std::cout << "\nUsage: full_cms [OPTIONS] INPUT_FILE\n\n" << std::endl;
  for (int i = 0; options[i].name != NULL; i++) {
    printf("\t-%c  --%s\t%s\n", options[i].val, options[i].name, options[i].has_arg ? options[i].name : "");
  }
  std::cout << "\n " << std::setw(100) << std::setfill('=') << "" << std::setfill(' ') << std::endl;
}

void GetInputArguments(int argc, char **argv)
{
  // process arguments
  if (argc == 1) {
    Help();
    exit(0);
  }
  while (true) {
    int c, optidx = 0;
    c = getopt_long(argc, argv, "pm:f:", options, &optidx);
    if (c == -1) break;
    //
    switch (c) {
    case 0:
      c = options[optidx].val;
    case 'p':
      parIsPerformance = true;
      break;
    case 'm':
      parMacroFileName = optarg;
      break;
    case 'f':
      parPhysListName = optarg;
      break;
    default:
      Help();
      errx(1, "unknown option %c", c);
    }
  }
  // check if mandatory Geant4 macro file was provided
  if (parMacroFileName == "") {
    G4cout << "  *** ERROR : Geant4 macro file is required. " << G4endl;
    Help();
    exit(-1);
  }
}
