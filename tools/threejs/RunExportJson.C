#include "TMath.h"
#include "TControlBar.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TVirtualPad.h"
#include "TCanvas.h"
#include "TVirtualGeoPainter.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TView.h"
#include "TPaveText.h"
#include "TGeoBBox.h"
#include "TGeoPara.h"
#include "TGeoTube.h"
#include "TGeoCone.h"
#include "TGeoEltu.h"
#include "TGeoSphere.h"
#include "TGeoTorus.h"
#include "TGeoTrd1.h"
#include "TGeoTrd2.h"
#include "TGeoParaboloid.h"
#include "TGeoHype.h"
#include "TGeoPcon.h"
#include "TGeoPgon.h"
#include "TGeoArb8.h"
#include "TGeoXtru.h"
#include "TGeoCompositeShape.h"
#include "TGeoPhysicalNode.h"
#include <algorithm>
#include <queue>  // std::queue
#include <math.h> /* sqrt */
#include <TGDMLWrite.h>

void RunExportJson()
{
  // Load some root geometry
  TGeoVolume *top   = gGeoManager->GetTopVolume();
  Int_t nTotalNodes = gGeoManager->GetNNodes();
  TGeoIterator iter(top);
  TGeoNode *current;
  TGeoVolume *vol;
  TGeoShape *shape;
  TBuffer3D *buffer_vol;
  TString path;
  TGeoVolume *CurrentVolume;
  TString guid;
  ofstream JsonFile;
  TString JsonFileName = "geo-box.json";
  TGeoVolume *MapVol[10000];
  TString MapGuidVol[10000];
  TGeoMaterial *MapMat[10000];
  TString MapGuidMat[10000];
  TString guidVol, guidMat;
  int iVol;

  // current=top;

  cout << "\n\n--> Starting elaboration of existing Geometry and creation of JSON file: " << JsonFileName << "\n";
  cout << "-->Top volume: " << top->GetName() << " Nodes: " << nTotalNodes << std::endl;

  JsonFile.open(JsonFileName);

  JsonFile << "{"
           << "\n"
           << "\"metadata\": { \n"
           << "\"version\": 4.3,"
           << "\"type\": \"Object\",\n"
           << "\"generator\": \"ObjectExporter\"},\n"
           << "\"geometries\": [\n";

  TIter next(gGeoManager->GetListOfVolumes());
  int kVolume = 0;
  cout << "Defining Volumes\n";
  while ((CurrentVolume = (TGeoVolume *)next())) {

    if (kVolume > 0) JsonFile << ",\n";
    shape         = CurrentVolume->GetShape();
    TString chvol = CheckVolume(shape);

    // cout <<"\n "<<kVolume<<"*****CHECK VOLUME: "<<"Class: "<<chvol<<"-->"<<CurrentVolume->GetName();
    // cout <<".";

    if (chvol != "") {
      if (strcmp(chvol, "TGeoCompositeShape") == 0) {
        cout << "\n ===Start: " << chvol << "\n";
        TGeoShape *leftS = GetLshape((TGeoCompositeShape *)shape);
        guid             = doGuidStuff(GuidGenerator());
        CreateJsonVol(leftS, JsonFile, guid);
        MapVol[kVolume]     = CurrentVolume;
        MapGuidVol[kVolume] = guid;
        kVolume++;
        JsonFile << ","
                 << "\n";
        TGeoShape *rightS = GetRshape((TGeoCompositeShape *)shape);
        guid              = doGuidStuff(GuidGenerator());
        CreateJsonVol(rightS, JsonFile, guid);
        MapVol[kVolume]     = CurrentVolume;
        MapGuidVol[kVolume] = guid;
        kVolume++;
        cout << "\n ===End: " << chvol << "\n";

      } else {
        guid = doGuidStuff(GuidGenerator());
        CreateJsonVol(shape, JsonFile, guid);
        MapVol[kVolume]     = CurrentVolume;
        MapGuidVol[kVolume] = guid;
        kVolume++;
      }
    }
  }
  JsonFile << "],"
           << "\n"; // End definition of Geometries
  JsonFile << "\"materials\": ["
           << "\n"; // Start definition of Materials

  // go through materials  - iterator and object declaration
  TIter nextmat(gGeoManager->GetListOfMaterials());
  TGeoMaterial *lmaterial;
  int kMaterial = 0;
  while ((lmaterial = (TGeoMaterial *)nextmat())) {
    if (kMaterial > 0) JsonFile << ",\n";
    // generate uniq name
    TString lname         = lmaterial->GetName();
    guid                  = doGuidStuff(GuidGenerator());
    MapMat[kMaterial]     = lmaterial;
    MapGuidMat[kMaterial] = guid;

    cout << "\n kMaterial :" << kMaterial << " " << lmaterial->GetName() << "\n";

    ExportMaterials(guid, lmaterial->GetName(), lmaterial->GetA(), lmaterial->GetZ(), JsonFile);

    kMaterial++;
  }
  // TOBE DELETED BEGIN

  // JsonFile << "}\n}";

  // JsonFile.close(); // close the file handle
  // std::cout << "\n\n--> JSON FILE CREATED "<< std::endl<< std::endl;
  // return;
  // TOBE DELETED END

  JsonFile << "],"
           << "\n"; // End definition of Materials

  JsonFile << "\"object\": {"
           << "\n"; // Start definition of scene
  guid = doGuidStuff(GuidGenerator());
  JsonFile << "\"uuid\": \"" << guid << "\",\n";
  JsonFile << "\"type\": \"Scene\",\n";
  JsonFile << "\"name\": \"Scene\",\n";
  JsonFile << "\"matrix\": [1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]";

  cout << "\n\n ==================Nodes====================\n\n";

  //  --Loop on the volumes of the tree structure
  int kLevel = 0;
  int kNode  = 0;
  while ((current = iter.Next())) {
    iter.GetPath(path);
    vol   = current->GetVolume();
    shape = vol->GetShape();

    TString chvol = CheckVolume(shape);
    // cout <<"\n *****CHECK Node class:--> "<< chvol<<", Vol. Name:-->"<<vol->GetName()<<"\n";
    if (chvol != "") {

      //   find the guid of the current volume
      guidVol = "";
      for (int i = 0; i < kVolume; i++) {
        if (vol == MapVol[i]) {
          iVol    = i;
          guidVol = MapGuidVol[i];
          // cout << "--"<<current->GetName()<<": ";
          break;
        }
      }
      if (guidVol == "") {
        cout << "ERROR **** UNDEFINED volume" << current->GetName() << "  STOP\n";
        break;
      }
      //   find the guid of the material of the current volume
      guidMat = "";
      for (int i = 0; i < kMaterial; i++) {
        if (vol->GetMaterial() == MapMat[i]) {
          // cout << vol->GetMaterial()->GetName();
          guidMat = MapGuidMat[i];
          break;
        }
      }
      if (guidMat == "") {
        guidMat = MapGuidMat[0];
        cout << "ERROR **** UNDEFINED material for:" << current->GetName() << "  ASSIGNED index 0\n";
      }

      if (iter.GetLevel() > kLevel) {
        //   If the level is greater that current level, then we define a new child
        JsonFile << ",\n";
        TabLine(iter.GetLevel(), JsonFile);
        JsonFile << "\"children\": [\n";
        TabLine(iter.GetLevel(), JsonFile);
        JsonFile << "{\n";
      } else if (iter.GetLevel() < kLevel) {
        int nLev = kLevel - iter.GetLevel();
        JsonFile << "\n";
        for (int i = 0; i < nLev; i++) {
          TabLine(kLevel - i, JsonFile);
          JsonFile << "}]\n";
        }
        TabLine(iter.GetLevel(), JsonFile);
        JsonFile << "},\n";
        TabLine(iter.GetLevel(), JsonFile);
        JsonFile << "{\n";

      } else {
        JsonFile << "\n";
        TabLine(iter.GetLevel(), JsonFile);
        JsonFile << "},\n";
        TabLine(iter.GetLevel(), JsonFile);
        JsonFile << "{\n";
      }

      if (strcmp(chvol, "TGeoCompositeShape") == 0) {
        cout << "\n ===Start Node: " << chvol << "\n";
        guid = doGuidStuff(GuidGenerator());
        ExportLeftCompVolume(guid, guidVol, guidMat, (TGeoCompositeShape *)shape, iter.GetLevel(), JsonFile);
        TabLine(iter.GetLevel(), JsonFile);
        JsonFile << "},\n";
        TabLine(iter.GetLevel(), JsonFile);
        JsonFile << "{\n";
        guidVol = MapGuidVol[iVol + 1];
        guid    = doGuidStuff(GuidGenerator());
        ExportRightCompVolume(guid, guidVol, guidMat, (TGeoCompositeShape *)shape, iter.GetLevel(), JsonFile);

        cout << "\n ===End Node: " << chvol << "\n";

      } else {

        guid = doGuidStuff(GuidGenerator());
        ExportCurrentVolume(guid, guidVol, guidMat, current, iter.GetLevel(), JsonFile);

        kNode++;
      }
    }
    kLevel = iter.GetLevel();
  }

  JsonFile << "\n";
  for (int i = 0; i < kLevel; i++) {
    TabLine(kLevel - i, JsonFile);
    JsonFile << "}]\n";
  }
  JsonFile << "}\n}";

  JsonFile.close(); // close the file handle

  std::cout << "\n\n--> JSON FILE: " << JsonFileName << " CREATED " << std::endl << std::endl;
}

//============== JSON GEOMETRY ======== End
