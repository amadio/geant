#ifndef GXRunManager_H
#define GXRunManager_H 1

#include "GXVCoprocessorManager.h"

#include <fstream>

using namespace std;

class GXRunManager 
{
public:
 
  static GXRunManager* Instance();
  GXRunManager();
  ~GXRunManager();

  void ConstructGeometry(GXVGeometry *geom);

  //only one coprocessor manager for now
  void RegisterCoprocessorManager(GXVCoprocessorManager *manager);

  void Initialization();
  void InitializeCoprocessorManager();

  GXVGeometry* GetGeometry();
  GXFieldMap** GetMagneticFieldMap();
  GXPhysicsTable* GetPhysicsTable();
  GXPhysics2DVector* GetSBData();

  inline GXVCoprocessorManager* GetCoprocessorManager() 
  { return coprocessorManager; };


private:
  
  void ReadMagneticFieldMap(const char *fieldMapFile); 
  void PreparePhysicsTable(); 
  void PreparePhysics2DVector(); 

  void readTable(GXPhysicsTable* table, const char* fname);
  void readTableAndSetSpline(GXPhysicsTable* table, const char* fname);
  bool RetrieveSeltzerBergerData(std::ifstream& in, GXPhysics2DVector *vec2D);

private:
  static GXRunManager* theInstance;

  GXVGeometry *geom;
  GXFieldMap** fieldMap;
  GXPhysicsTable* physicsTable;
  GXPhysics2DVector* sbData;

  GXTaskData* taskData;

  GXVCoprocessorManager* coprocessorManager;
};

#endif
