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

//============== GUID =================begin
#include "/Users/goulas/Desktop/geom/guid.h"
/*
The MIT License (MIT)

Copyright (c) 2014 Graeme Hill (http://graemehill.ca)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#include "./guid.h"

#ifdef GUID_LIBUUID
#include <uuid/uuid.h>
#endif

#include <CoreFoundation/CFUUID.h>

using namespace std;

// overload << so that it's easy to convert to a string
ostream &operator<<(ostream &s, const Guid &guid)
{
  return s << hex << setfill('0') << setw(2) << (int)guid._bytes[0] << setw(2) << (int)guid._bytes[1] << setw(2)
           << (int)guid._bytes[2] << setw(2) << (int)guid._bytes[3] << "-" << setw(2) << (int)guid._bytes[4] << setw(2)
           << (int)guid._bytes[5] << "-" << setw(2) << (int)guid._bytes[6] << setw(2) << (int)guid._bytes[7] << "-"
           << setw(2) << (int)guid._bytes[8] << setw(2) << (int)guid._bytes[9] << "-" << setw(2) << (int)guid._bytes[10]
           << setw(2) << (int)guid._bytes[11] << setw(2) << (int)guid._bytes[12] << setw(2) << (int)guid._bytes[13]
           << setw(2) << (int)guid._bytes[14] << setw(2) << (int)guid._bytes[15];
}

// create a guid from vector of bytes
Guid::Guid(const vector<unsigned char> &bytes)
{
  _bytes = bytes;
}

// create a guid from array of bytes
Guid::Guid(const unsigned char *bytes)
{
  _bytes.assign(bytes, bytes + 16);
}

// converts a single hex char to a number (0 - 15)
unsigned char hexDigitToChar(char ch)
{
  if (ch > 47 && ch < 58) return ch - 48;

  if (ch > 96 && ch < 103) return ch - 87;

  if (ch > 64 && ch < 71) return ch - 55;

  return 0;
}

// converts the two hexadecimal characters to an unsigned char (a byte)
unsigned char hexPairToChar(char a, char b)
{
  return hexDigitToChar(a) * 16 + hexDigitToChar(b);
}

// create a guid from string
Guid::Guid(const string &fromString)
{
  _bytes.clear();

  char charOne, charTwo;
  bool lookingForFirstChar = true;

  for (const char &ch : fromString) {
    if (ch == '-') continue;

    if (lookingForFirstChar) {
      charOne             = ch;
      lookingForFirstChar = false;
    } else {
      charTwo   = ch;
      auto byte = hexPairToChar(charOne, charTwo);
      _bytes.push_back(byte);
      lookingForFirstChar = true;
    }
  }
}

// create empty guid
Guid::Guid()
{
  _bytes = vector<unsigned char>(16, 0);
}

// copy constructor
Guid::Guid(const Guid &other)
{
  _bytes = other._bytes;
}

// overload assignment operator
Guid &Guid::operator=(const Guid &other)
{
  _bytes = other._bytes;
  return *this;
}

// overload equality operator
bool Guid::operator==(const Guid &other) const
{
  return _bytes == other._bytes;
}

// overload inequality operator
bool Guid::operator!=(const Guid &other) const
{
  return !((*this) == other);
}

Guid GuidGenerator::newGuid()
{
  auto newId = CFUUIDCreate(NULL);
  auto bytes = CFUUIDGetUUIDBytes(newId);
  CFRelease(newId);

  const unsigned char byteArray[16] = {
      bytes.byte0, bytes.byte1, bytes.byte2,  bytes.byte3,  bytes.byte4,  bytes.byte5,  bytes.byte6,  bytes.byte7,
      bytes.byte8, bytes.byte9, bytes.byte10, bytes.byte11, bytes.byte12, bytes.byte13, bytes.byte14, bytes.byte15};
  return byteArray;
}

TString UpperCase(TString str)
{
  for (int i = 0; i < strlen(str); i++)
    str[i]   = toupper(str[i]);
  return str;
}

TString doGuidStuff(GuidGenerator generator)
{
  auto myGuid = generator.newGuid();
  std::stringstream stream;
  stream << myGuid;
  auto guidString = stream.str();

  return UpperCase((TString)guidString);
}

//============== GUID ================= End

//============== JSON GEOMETRY ======== Begin
//______________________________________________________________________________

void CalculateSurfaceNormal(Double_t p1[3], Double_t p2[3], Double_t p3[3], Double_t *Vnormal)
{
  //
  // Here we compute two vectors coplanar to a face. The faces are triangles defined by points
  // p1,p2,p3, then two possible vectors are: U = p2 - p1 and V = p3 - p1
  // With the two vectors, U & V we compute the cross product between them to find a perpendicular
  // vector to the face: Vn = U x V
  // Vn(x) = Uy * Vz - Uz * Vy
  // Vn(y) = Uz * Vx - Ux * Vz
  // Vn(z) = Ux * Vy - Uy * Vx
  //
  // The length of the vector: Vmod=SQRT(Vnx**2+Vny**2+Vnz**2)
  //
  //
  Double_t U[3], V[3], Vn[3], Vmod;

  U[0] = p2[0] - p1[0];
  U[1] = p2[1] - p1[1];
  U[2] = p2[2] - p1[2];

  V[0] = p3[0] - p1[0];
  V[1] = p3[1] - p1[1];
  V[2] = p3[2] - p1[2];

  Vn[0] = U[1] * V[2] - U[2] * V[1];
  Vn[1] = U[2] * V[0] - U[0] * V[2];
  Vn[2] = U[0] * V[1] - U[1] * V[0];

  Vmod = sqrt(pow(Vn[0], 2.0) + pow(Vn[1], 2.0) + pow(Vn[2], 2.0));

  if (Vmod == 0) {
    cout << "\n-->Error:====================================================================\n";
    for (int i = 0; i < 3; i++) {
      cout << "-->P: " << p1[i] << " " << p2[i] << " " << p3[i] << " U: " << U[i] << " V: " << V[i] << " Vn: " << Vn[i]
           << "  ||  ";
    }
  }

  Vnormal[0] = Vn[0] / Vmod;
  Vnormal[1] = Vn[1] / Vmod;
  Vnormal[2] = Vn[2] / Vmod;
  // cout << "Vmode:" << Vmod <<" Normal:" << Vnormal <<"\n";
}

void CreateJsonVol(TGeoVolume *vol, TGeoShape *shape, ofstream &JsonFile, TString guid)
{
  // Load some root geometry
  TGeoNode *current;
  TBuffer3D *buffer_vol;
  Nfaces = 0;
  queue<int *> Faces;
  int CPh[200000][3];
  int *FaceTriangle;
  queue<Double_t *> Fnormals;
  Double_t V1[3], V2[3], V3[3];
  Double_t *Vnormal;
  bool found;

  //--Consider each polygon of min 3 (triangles) and max 4 (quads) vertices.
  //--vrtx containes the index to vertexes for each segment of the polugon (max 2*4=8)
  //--vquad contains [continous] vertexes of the polygon (max=4)
  int vrtx[8], vquad[4], indxf;

  //     --Call the triangulation routine
  buffer_vol = shape->MakeBuffer3D();

  cout << "\n-->No of Vertices: " << buffer_vol->NbPnts() << "\n";
  cout << "-->No of Segments: " << buffer_vol->NbSegs() << "\n";
  cout << "-->No of Polygons: " << buffer_vol->NbPols() << "\n";
  cout << "\n-->Writing " << buffer_vol->NbPnts() << " Vertices to JSON file ... ";

  // Start GEOMETRY definition in JSON with "{".
  JsonFile << "{\n";
  JsonFile << "\t\"uuid\": \"";
  // Generate a GUID for the specific volume for JSON file
  JsonFile << guid << "\",\n";
  JsonFile << "\t\"type\": \"Geometry\",\n";
  JsonFile << "\t\"data\":{\n";
  JsonFile << "\t\"vertices\": [";
  // Write all vertexes in JSON file
  int Kvertex = 0;
  for (int i = 0; i < buffer_vol->NbPnts() * 3; i = i + 3) {
    JsonFile << buffer_vol->fPnts[i] << ",  " << buffer_vol->fPnts[i + 1] << ",  " << buffer_vol->fPnts[i + 2];
    if (i < buffer_vol->NbPnts() * 3 - 3) {
      JsonFile << ","; // JsonFile <<",\n";
    } else {
      // JsonFile <<"\n";
    }
    Kvertex = Kvertex + 1;
  }
  JsonFile << "],"
           << "\n";
  cout << "Done vertex";

  IndexSegPolygon = 1; // define only triangles
  Nfaces          = 0;
  cout << "\n-->Generating triangles and normals ...";

  //      Loop through all the polygons that have been defined for this volume
  //      Each polygon will be converted to triangles faces saved in the que <Faces>
  //      and for each triangle we calculate the normal stored also in the queue <Fnormals>

  for (int Kpolygon = 0; Kpolygon < buffer_vol->NbPols(); Kpolygon++) {
    int iv = 0;

    //          Loop through all the segments of the current polygon and get the Vertex indexes
    //          expected max 4 segments and consequently 8 vertex indexes
    //          fSegs:  c0, p0, q0, c1, p1, q1, ..... ..... ....
    //          fPols;  c0, n0, s0, s1, ... sn, c1, n1, s0, ... sn

    for (int k = 1; k < buffer_vol->fPols[IndexSegPolygon] + 1; k++) {
      vrtx[iv] = buffer_vol->fSegs[3 * buffer_vol->fPols[IndexSegPolygon + k] + 1];
      iv       = iv + 1;
      vrtx[iv] = buffer_vol->fSegs[3 * buffer_vol->fPols[IndexSegPolygon + k] + 2];
      iv       = iv + 1;
    }

    //          Segments are defined in the correct order for the definition of the normals
    //          but for each segment the sequence of vertexes may be inverted.
    //          check continuity of vertexes in the consequtive segments
    //          If Seg1 made of P0-P1 and Seg2 made of P0-P2 then Seg1 should be defined as P1-P0

    if (vrtx[1] != vrtx[2]) {
      if (vrtx[0] == vrtx[2]) {
        int itemp = vrtx[0];
        vrtx[0]   = vrtx[1];
        vrtx[1]   = itemp;

      } else if (vrtx[0] == vrtx[3]) {
        int itemp = vrtx[0];
        vrtx[0]   = vrtx[1];
        vrtx[1]   = itemp;

        itemp   = vrtx[2];
        vrtx[2] = vrtx[3];
        vrtx[3] = itemp;

      } else if (vrtx[1] == vrtx[3]) {
        int itemp = vrtx[2];
        vrtx[2]   = vrtx[3];
        vrtx[3]   = itemp;
      }
    }

    //          ---
    //          Out of the (max 8) vertexes definining the (max 4) segments,
    //          get the (max 4) vertexes defining the polygon (in quad)
    //          Check that segments and vertexes are in the correct order
    //          Starting from the first segment (vrtx[0]-vrtx[1]),
    //          the second segment must have the vertex vrtx[1] and this should be the first in the order
    //          ---

    vquad[0] = vrtx[0];
    vquad[1] = vrtx[1];
    vrtx[0]  = -1;
    vrtx[1]  = -1;

    int kvert = 2;
    while (kvert < buffer_vol->fPols[IndexSegPolygon]) {
      for (int ks = 1; ks < iv; ks++) {
        if (vrtx[ks * 2 + 1] == vquad[kvert - 1]) {
          vquad[kvert]     = vrtx[ks * 2];
          vrtx[ks * 2]     = -1;
          vrtx[ks * 2 + 1] = -1;
          kvert            = kvert + 1;
          break;
        } else if (vrtx[ks * 2] == vquad[kvert - 1]) {
          vquad[kvert]     = vrtx[ks * 2 + 1];
          vrtx[ks * 2]     = -1;
          vrtx[ks * 2 + 1] = -1;
          kvert            = kvert + 1;
          break;
        }
      }
    }

    //     Using the vertex indexes get the x,yz coordinates of each vertex of the 1st triangle

    for (int i = 0; i < 3; i++) {
      V1[i] = buffer_vol->fPnts[3 * vquad[2] + i];
      V2[i] = buffer_vol->fPnts[3 * vquad[1] + i];
      V3[i] = buffer_vol->fPnts[3 * vquad[0] + i];
    }

    if ((V1[0] == V2[0] && V1[1] == V2[1] && V1[2] == V2[2]) || (V1[0] == V3[0] && V1[1] == V3[1] && V1[2] == V3[2]) ||
        (V2[0] == V3[0] && V2[1] == V3[1] && V2[2] == V3[2])) {

      cout << "\nSkiping triangle face: <" << Nfaces << "> [" << V1[0] << "," << V1[1] << " ," << V1[2] << "]";

    } else {

      //            Get the first triangle of the defined polygon

      FaceTriangle    = new int[3];
      FaceTriangle[0] = vquad[2];
      FaceTriangle[1] = vquad[1];
      FaceTriangle[2] = vquad[0];

      Vnormal = new Double_t[3];
      CalculateSurfaceNormal(V1, V2, V3, Vnormal);

      Fnormals.push(Vnormal);

      Faces.push(FaceTriangle);
      Nfaces = Nfaces + 1;
    }

    //          If the polygon is a quad, get also the second triangle of the defined polygon
    //
    //         3 _____ 2
    //          |    /|
    //          | B / |
    //          |  /  |
    //          | / A |
    //          |/____|
    //         0       1
    //
    //         First trianle A: 0,1,2 second triangle B: 2,3,0

    if (buffer_vol->fPols[IndexSegPolygon] == 4) {

      //            Using the vertex indexes get the x,yz coordinates of each vertex of the 2nd triangle

      for (int i = 0; i < 3; i++) {
        V1[i] = buffer_vol->fPnts[3 * vquad[0] + i];
        V2[i] = buffer_vol->fPnts[3 * vquad[3] + i];
        V3[i] = buffer_vol->fPnts[3 * vquad[2] + i];
      }

      if ((V1[0] == V2[0] && V1[1] == V2[1] && V1[2] == V2[2]) ||
          (V1[0] == V3[0] && V1[1] == V3[1] && V1[2] == V3[2]) ||
          (V2[0] == V3[0] && V2[1] == V3[1] && V2[2] == V3[2])) {

        cout << "\nSkiping triangle face: <" << Nfaces << "> [" << V1[0] << "," << V1[1] << " ," << V1[2] << "]";

      } else {

        FaceTriangle    = new int[3];
        FaceTriangle[0] = vquad[0];
        FaceTriangle[1] = vquad[3];
        FaceTriangle[2] = vquad[2];

        Vnormal = new Double_t[3];
        CalculateSurfaceNormal(V1, V2, V3, Vnormal);

        Fnormals.push(Vnormal);

        Faces.push(FaceTriangle);

        Nfaces = Nfaces + 1;
      }
    }

    IndexSegPolygon = IndexSegPolygon + buffer_vol->fPols[IndexSegPolygon] + 2;
  } // end loop through polygons

  //      print the Normals (one normal per triangle)

  JsonFile << "\t\"normals\": [";
  cout << "\n-->Generated " << Nfaces << " Normals ... Writing normals to file ... ";
  int Knormal = 0;
  while (!Fnormals.empty()) {
    Double_t *Vnorm = Fnormals.front();
    JsonFile << Vnorm[0] << "," << Vnorm[1] << "," << Vnorm[2];
    Knormal = Knormal + 1;

    if (Nfaces != Knormal) {
      JsonFile << ",";
    }
    // JsonFile <<" \n";
    delete[] Vnorm;
    Fnormals.pop();
  }
  JsonFile << " ],\n";

  // Start writing faces definition in Json file

  JsonFile << "\t\"faces\": [";

  cout << "Done Normals\n";
  cout << "-->Creating triangle faces ... ";

  //      print the triangles / faces definition
  int knorm  = 0;
  int Kfaces = 0;
  while (!Faces.empty()) {
    int *point = Faces.front();
    JsonFile << "16," << point[0] << "," << point[1] << "," << point[2] << "," << knorm << "  ";
    CPh[Kfaces][0] = point[0];
    CPh[Kfaces][1] = point[1];
    CPh[Kfaces][2] = point[2];
    Kfaces         = Kfaces + 1;
    knorm          = knorm + 1;

    if (Nfaces != Kfaces) {
      JsonFile << ",";
    }
    // JsonFile <<"\n";
    delete[] point;
    Faces.pop();
  }
  JsonFile << "]";
  // Terminate the DATA definition for this volume with "}".
  JsonFile << "}\n";

  // Terminate the GEOMETRY definition for this volume with "}".
  JsonFile << "}";

  cout << Nfaces << " triangles created\n";
  cout << "-->Checking vertex order for each triangle ... ";

  //      Check the def. of the triangles if the order of vertixes is correct

  Kfaces = 0;
  while (Kfaces < Nfaces) {
    int Kf = 0;
    found  = false;
    while (Kf < Nfaces) {
      if (Kf != Kfaces) {
        if (((CPh[Kfaces][0] == CPh[Kf][1]) && (CPh[Kfaces][1] == CPh[Kf][0])) ||
            ((CPh[Kfaces][0] == CPh[Kf][2]) && (CPh[Kfaces][1] == CPh[Kf][1])) ||
            ((CPh[Kfaces][0] == CPh[Kf][0]) && (CPh[Kfaces][1] == CPh[Kf][2])) ||

            ((CPh[Kfaces][1] == CPh[Kf][1]) && (CPh[Kfaces][2] == CPh[Kf][0])) ||
            ((CPh[Kfaces][1] == CPh[Kf][2]) && (CPh[Kfaces][2] == CPh[Kf][1])) ||
            ((CPh[Kfaces][1] == CPh[Kf][0]) && (CPh[Kfaces][2] == CPh[Kf][2])) ||

            ((CPh[Kfaces][2] == CPh[Kf][1]) && (CPh[Kfaces][0] == CPh[Kf][0])) ||
            ((CPh[Kfaces][2] == CPh[Kf][2]) && (CPh[Kfaces][0] == CPh[Kf][1])) ||
            ((CPh[Kfaces][2] == CPh[Kf][0]) && (CPh[Kfaces][0] == CPh[Kf][2]))) {
          found = true;
          cout << ".";
          break;
        }
      }
      Kf++;
    }
    if (!found) {
      cout << "<Error on face:" << Kfaces << ">";
    }
    Kfaces++;
  }
  cout << " Done Checking faces\n";
}

//============== JSON Materials definition ======== Begin
void ExportMaterials(TString guid, TString MatName, Double_t nA, Double_t nZ, ofstream &JsonFile)
{
  JsonFile << "{\n\t\"uuid\": \"" << guid << "\",\n";
  JsonFile << "\t\"type\": \"MeshLambertMaterial\",\n";
  JsonFile << "\t\"name\":\"" << MatName << "\",\n";
  JsonFile << "\t\"color\":16777215,\n";
  JsonFile << "\t\"emissive\": 16723858,\n";
  if (MatName == "Vacuum" || (nA == 0 && nZ == 0)) {
    JsonFile << "\t\"opacity\": 0.0,\n";
    JsonFile << "\t\"transparent\": true\n }";
  } else {
    JsonFile << "\t\"opacity\": 1.0,\n";
    JsonFile << "\t\"transparent\": false\n }";
  }
}

//============== JSON Materials definition ======== End

void TabLine(int ntab, ofstream &JsonFile)
{
  for (i = 0; i < ntab; i++) {
    JsonFile << "\t";
  }
}

//============== JSON Nodes definition ======== Begin

TString GetHexRGB(int r, int g, int b)
{
  char hexcol[16];

  snprintf(hexcol, sizeof hexcol, "%02x%02x%02x", r, g, b);
  return hexcol;
}

void scan(TString guid, TString guidVol, TString guidMat, TString VolNane, int kLevel, ofstream &JsonFile)
{
  TabLine(kLevel + 1, JsonFile);
  JsonFile << "\"uuid\": \"" << guid << "\",\n";
  TabLine(kLevel + 1, JsonFile);
  JsonFile << "\"type\": \"Mesh\",\n";
  TabLine(kLevel + 1, JsonFile);
  JsonFile << "\"name\": \"" << VolNane << "\",\n";
  TabLine(kLevel + 1, JsonFile);
  JsonFile << "\"geometry\": \"" << guidVol << "\",\n";
  TabLine(kLevel + 1, JsonFile);
  JsonFile << "\"material\": \"" << guidMat << "\",\n";
  TabLine(kLevel + 1, JsonFile);
  JsonFile << "\"matrix\": [1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]";
}
//============== JSON Nodes definition ======== End

void ExportJson()
{
  // Load some root geometry
  TGeoVolume *top = gGeoManager->GetTopVolume();
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

  // FILE* file = fopen(JsonFileName, "wb");
  // open a file for writing the son file format
  cout << "\n\n--> Starting elaboration of existing Geometry and cration of JSON file: " << JsonFileName << " \n\n";
  cout << "--> JSON data in file: " << JsonFileName << ">\n";

  cout << "-->Top volume: " << top->GetName() << std::endl;

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
  while ((CurrentVolume = (TGeoVolume *)next())) {
    if (kVolume > 0) JsonFile << ",\n";
    printf("Current volume is: %s\n", CurrentVolume->GetName());
    guid  = doGuidStuff(GuidGenerator());
    shape = CurrentVolume->GetShape();
    CreateJsonVol(CurrentVolume, shape, JsonFile, guid);
    MapVol[kVolume]     = CurrentVolume;
    MapGuidVol[kVolume] = guid;
    cout << MapVol[kVolume] << "  --  " << MapGuidVol[kVolume];
    kVolume++;
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

    ExportMaterials(guid, lmaterial->GetName(), lmaterial->GetA(), lmaterial->GetZ(), JsonFile);
    kMaterial++;
  }
  JsonFile << "],"
           << "\n"; // End definition of Materials

  JsonFile << "\"object\": {"
           << "\n"; // Start definition of scene
  guid = doGuidStuff(GuidGenerator());
  JsonFile << "\"uuid\": \"" << guid << "\",\n";
  JsonFile << "\"type\": \"Scene\",\n";
  JsonFile << "\"name\": \"Scene\",\n";
  JsonFile << "\"matrix\": [1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]";

  //  --Loop on the volumes of the tree structure
  int kLevel = 0;
  int kNode  = 0;
  while ((current = iter.Next())) {
    iter.GetPath(path);
    vol   = current->GetVolume();
    shape = vol->GetShape();
    //   find the guid of the current volume
    guidVol = "";
    for (int i = 0; i < kVolume; i++) {
      if (vol == MapVol[i]) {
        guidVol = MapGuidVol[i];
        cout << "--" << current->GetName() << ": ";
        break;
      }
    }
    if (guidVol == "") {
      cout << "ERROR **** UNDEFINED volume" << current->GetName() << "  STOP\n";
      // break;
    }
    //   find the guid of the material of the current volume
    guidMat = "";
    for (int i = 0; i < kMaterial; i++) {
      if (vol->GetMaterial() == MapMat[i]) {
        cout << vol->GetMaterial()->GetName();
        guidMat = MapGuidMat[i];
        break;
      }
    }
    if (guidMat == "") {
      cout << "ERROR **** UNDEFINED material for:" << current->GetName() << "  STOP\n";
      // break;
    }
    if (iter.GetLevel() > kLevel) {
      //   If the level is greater that current level, then we define a new child
      JsonFile << ",\n";
      TabLine(iter.GetLevel(), JsonFile);
      JsonFile << "\"children\": [\n";
      TabLine(iter.GetLevel(), JsonFile);
      JsonFile << "{\n";
      //      TabLine(iter.GetLevel()+1,JsonFile);
      //      JsonFile << current->GetName()<<"  level:"<<iter.GetLevel()<<" -> "<< kLevel<<"\n";
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
      //      TabLine(iter.GetLevel()+1,JsonFile);
      //      JsonFile << current->GetName()<<"  level:"<<iter.GetLevel()<<" -> "<< kLevel<<"\n";

    } else {
      JsonFile << "\n";
      TabLine(iter.GetLevel(), JsonFile);
      JsonFile << "},\n";
      TabLine(iter.GetLevel(), JsonFile);
      JsonFile << "{\n";
      //      TabLine(iter.GetLevel()+1,JsonFile);
      //      JsonFile << current->GetName()<<"  level:"<<iter.GetLevel()<<" -> "<< kLevel<<"\n";
    }
    kLevel = iter.GetLevel();

    guid = doGuidStuff(GuidGenerator());

    scan(guid, guidVol, guidMat, current->GetName(), kLevel, JsonFile);

    // cout << "\n-->Current path: " << path << " Name: " << current->GetName() <<" Level: "<< iter.GetLevel()<< " <"<<
    // std::hex <<GetHexRGB(0xff, 0xff, 0xff) <<">"<<std::endl;
    kNode++;
  }
  JsonFile << "\n";
  for (int i = 0; i < kLevel; i++) {
    TabLine(kLevel - i, JsonFile);
    JsonFile << "}]\n";
  }
  JsonFile << "}\n}";

  JsonFile.close(); // close the file handle
  std::cout << "\n\n--> JSON FILE CREATED " << std::endl << std::endl;
  // scan();

  JsonFile << "]\n";
}

//============== JSON GEOMETRY ======== End
