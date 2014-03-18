#include <math.h>
#include <iomanip>
#include <iostream>

#include "GXRunManager.h"
#include "GXUserGeometry.h"
#include "GXFieldMapData.h"
#include "GXPhysicsTableType.h"

using namespace std;

GXRunManager* GXRunManager::theInstance = 0;

GXRunManager* GXRunManager::Instance() {
  if(theInstance == 0) theInstance = new GXRunManager();
  return theInstance;
}

GXRunManager::GXRunManager()
{
  //geometry
  geom = NULL;

  //field map
  fieldMap = (GXFieldMap **) malloc (nbinZ*sizeof(GXFieldMap *));
  for (int j = 0 ; j < nbinZ ; j++) {
    fieldMap[j] = (GXFieldMap *) malloc (nbinR*sizeof(GXFieldMap));
  }

 //physics tables
  physicsTable = 
    (GXPhysicsTable*) malloc (kNumberPhysicsTable*sizeof(GXPhysicsTable));
 
  //DB data
  sbData = (GXPhysics2DVector*) malloc(maxElements*sizeof(GXPhysics2DVector));
}

GXRunManager::~GXRunManager()
{
  delete geom;
  free(fieldMap);
  free(physicsTable);
  free(sbData);
}

void GXRunManager::Initialization()
{ 
  //1. Construct geometry
  if(geom !=NULL)  geom->create();

  //2. Read magnetic field map
  const char* fieldMapFile = getenv("GP_BFIELD_MAP");
  fieldMapFile = (fieldMapFile) ? fieldMapFile : "data/cmsExp.mag.3_8T";
  //  GXFieldMap** fieldMap = ReadMagneticFieldMap(fieldMapFile);
  ReadMagneticFieldMap(fieldMapFile);

  //3. Prepare EM physics tables
  PreparePhysicsTable();
  PreparePhysics2DVector(); 

  //4. Initialize coprocessor manager
  InitializeCoprocessorManager();
}

void GXRunManager::ConstructGeometry(GXVGeometry *userDetector)
{ 
  geom = userDetector;
}

void GXRunManager::RegisterCoprocessorManager(GXVCoprocessorManager *manager)
{ 
  coprocessorManager = manager;
}

void GXRunManager::InitializeCoprocessorManager() 
{
  if (strcmp(coprocessorManager->GetName(),"GPU") == 0 ) {
    std::cout << "... Using GPU as a Coprocessor " << std::endl;
    coprocessorManager->Initialize();
    //    coprocessorManager->AllocateDeviceMemory(taskData);
    coprocessorManager->AllocateDeviceMemory(geom,fieldMap,
    					     physicsTable,sbData);
  }
}

//GXFieldMap** GXRunManager::ReadMagneticFieldMap(const char *fieldMapFile) 
void GXRunManager::ReadMagneticFieldMap(const char *fieldMapFile) 
{
  //@@@G4FWP: this method is a temporary solution and we should implement 
  //a generic interface for MagneticField

  // Read magnetic field map                                                  
  std::ifstream ifile(fieldMapFile, ios::in | ios::binary | ios::ate);

  if (ifile.is_open()) {

    //field map structure                                                       
    GXFieldMapData fd;

    ifstream::pos_type fsize = ifile.tellg();
    size_t dsize = sizeof(GXFieldMapData);

    long int ngrid = fsize/dsize;
    ifile.seekg (0, ios::beg);

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
  }
}

void GXRunManager::PreparePhysicsTable()
{
  //store physics tables

  char filename[256];
  for(int it = 0 ; it < kNumberPhysicsTable ; ++it) {
    sprintf(filename,"data/%s",GXPhysicsTableName[it]);
    readTableAndSetSpline(&physicsTable[it],filename);
  }
}

void GXRunManager::PreparePhysics2DVector()
{
  //store SeltzerBerger Data files

  char sbDataFile[256];
  for(G4int iZ = 0 ; iZ < maxElements; iZ++) {  
    sprintf(sbDataFile,"data/brem_SB/br%d",iZ+1);
    std::ifstream fin(sbDataFile);
    G4bool check = RetrieveSeltzerBergerData(fin, &sbData[iZ]);
    if(!check) {
      printf("Failed To open SeltzerBerger Data file for Z= %d\n",iZ+1);
    }
  }
}  

void GXRunManager::readTable(GXPhysicsTable* table, const char* fname) {
  std::ifstream fIn;
  fIn.open(fname,std::ios::in);
  if(!fIn){
    std::cout << fname << " was not found!!!" << std::endl;
    return;
  }
  fIn >> table->nPhysicsVector;
  for(int idx=0; idx<table->nPhysicsVector; idx++){
    int vType=0;
    fIn >> vType;
    table->physicsVectors[idx].type = vType;
    fIn >> table->physicsVectors[idx].edgeMin;
    fIn >> table->physicsVectors[idx].edgeMax;
    fIn >> table->physicsVectors[idx].numberOfNodes;
    int siz=0;
    fIn >> siz;
    for(int j=0; j<siz; j++){
      fIn >> table->physicsVectors[idx].binVector[j];
      fIn >> table->physicsVectors[idx].dataVector[j];
    }// j
    G4double theEmin = table->physicsVectors[idx].binVector[0];
    table->physicsVectors[idx].dBin = 
      log10(table->physicsVectors[idx].binVector[1]/theEmin);
    table->physicsVectors[idx].baseBin = 
      log10(theEmin)/table->physicsVectors[idx].dBin;
  }// idx
}

void GXRunManager::readTableAndSetSpline(GXPhysicsTable* table, 
					 const char* fname)
{
  bool useSpline = true;
  readTable(table,fname);
  G4int nv = table->nPhysicsVector;
  for(int j=0; j < nv; j++){
    table->physicsVectors[j].SetSpline(useSpline);
  }
}

bool GXRunManager::RetrieveSeltzerBergerData(std::ifstream& in, 
					     GXPhysics2DVector *vec2D)
{
  // binning
  G4int k;
  G4int dummyX; // 32 fixed up to Z = 92
  G4int dummyY; // 57 fixed up to Z = 92
  //  in >> k >> numberOfXNodes >> numberOfYNodes;
  in >> k >> dummyX >> dummyY;
  if (in.fail())  { return false; }

  // contents
  G4double valx, valy, val;
  for(size_t i = 0; i<numberOfXNodes; ++i) {
    in >> valx;
    if (in.fail())  { return false; }
    vec2D->PutX(i,valx);
   }
  for(size_t j = 0; j<numberOfYNodes; ++j) {
    in >> valy;
    if (in.fail())  { return false; }
    vec2D->PutY(j,valy);
   }
  for(size_t j = 0; j<numberOfYNodes; ++j) {
    for(size_t i = 0; i<numberOfXNodes; ++i) {
      in >> val;
      if (in.fail())  { return false; }
      vec2D->PutValue(i, j, val);
     }
  }
  in.close();
  return true;
}

GXVGeometry* GXRunManager::GetGeometry() {
  return geom;
}

GXFieldMap** GXRunManager::GetMagneticFieldMap()
{
  return fieldMap;
}

GXPhysicsTable* GXRunManager::GetPhysicsTable()
{
  return physicsTable;
}

GXPhysics2DVector* GXRunManager::GetSBData()
{
  return sbData;
}
