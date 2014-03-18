#include "cmsExpMagneticField.hh"
//#include "cmsExpTrackerBfield.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"

#include <sstream>
#include <iostream>
#include <fstream>

cmsExpMagneticField::cmsExpMagneticField()
  : flux(0.*tesla)
{
  //default 
  fieldType = cmsExp::kNull ;

  //prepare fieldMap array: fieldMap[cmsExp::nbinZ][cmsExp::nbinR];
  fieldMap = (G4TwoVector **) malloc (cmsExp::nbinZ*sizeof(G4TwoVector *));
  for (int j = 0 ; j < cmsExp::nbinZ ; j++) {
    fieldMap[j] = (G4TwoVector *) malloc (cmsExp::nbinR*sizeof(G4TwoVector));
  }

  const char* fieldMapFile = getenv("CMSEXP_BFIELD_MAP");
  fieldMapFile = (fieldMapFile) ? fieldMapFile : "cmsExp.mag.3_8T";

  ReadFieldMap(fieldMapFile);

  // theTrakcerBfield = new cmsExpTrackerBfield();

  // GL: taken from basic/B2/B2a example - all 3 lines necessary, otherwise there is no magnetic field
  // G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  // fieldMgr->SetDetectorField(this);
  // fieldMgr->CreateChordFinder(this);
}

cmsExpMagneticField::~cmsExpMagneticField()
{
  //delete theTrakcerBfield;
  //  free fieldMap;
}

void cmsExpMagneticField::GetFieldValue(const double point[3], 
					double *bField) const
{
  //default values
  bField[0] = 0.; 
  bField[1] = 0.;
  bField[2] = 0.;

  G4double rho = std::sqrt(point[0]*point[0]+point[1]*point[1]); 

  if(IsDefined(point[2],rho)) {

    switch(fieldType) {
      case cmsExp::kVolumebase :
	GetVolumeBaseBfield(point, bField);
	break;
      case cmsExp::kParametric :
	GetParametricBfield(point, bField, rho);
	break;
      case cmsExp::kUniform :
	GetUniformBfield(point, bField, rho);
	break;
      default :
	G4cout << "GetFieldType: Invalid BFieldType " << G4endl; 
    }   
  }

 }

void cmsExpMagneticField::SetFieldType(G4String newType) {

  if(newType == "volumebase") {
    fieldType = cmsExp::kVolumebase;
  }
  else if(newType == "parametric") {
    fieldType = cmsExp::kParametric;
  }
  else if(newType == "uniform") {
    fieldType = cmsExp::kUniform;
  }
  else {
    G4cout << "SetFieldType: Magnetic FieldType is not valid" << G4endl;
  }   
}

void cmsExpMagneticField::ReadFieldMap(const char* filename) {

  ifstream ifile(filename, ios::in | ios::binary | ios::ate);
  
  if (ifile.is_open()) {

    //field map structure
    cmsExp::cmsFieldMapData fd;

    ifstream::pos_type fsize = ifile.tellg();
    size_t dsize = sizeof(cmsExp::cmsFieldMapData);    

    long int ngrid = fsize/dsize;
    ifile.seekg (0, ios::beg);
    
    G4cout << "cmsExp ... Loading magnetic field map: " << filename << G4endl;

    for(G4int i = 0 ; i < ngrid ; i++) {
      ifile.read((char *)&fd, sizeof(cmsExp::cmsFieldMapData));
      
      //check validity of input data
      if(abs(fd.iz) > cmsExp::noffZ || fd.ir > cmsExp::nbinR) {
	G4Exception("cmsExpMagneticField::ReadFieldMap()", "BadIndex",
		    FatalException, "Out of magnetic field index!");
      }
      else {
	fieldMap[cmsExp::noffZ+fd.iz][fd.ir].set(fd.Bz,fd.Br);
      }
    }
    ifile.close();
  }
}  

bool cmsExpMagneticField::IsDefined(const double& z, const double& rho) const {
  return (std::fabs(z) < cmsExp::maxZ && rho < cmsExp::maxR);
}

void cmsExpMagneticField::GetVolumeBaseBfield(G4double const *point, 
					      G4double *bField)  const 
{
  //unit conversion: [mm] to [cm]
  G4double rho = 0.1*std::sqrt(point[0]*point[0]+point[1]*point[1]); 
  G4double   z = 0.1*point[2];

  //volume based magnetic field map - [0:900][-1600:1600] grid in [r][z]
  G4int     ir = int(rho);
  G4double  dr = rho - ir;

  // The previous implementation was:
  //    G4int iz = int(point[2]) + cmsExp::noffZ;
  // the main issue is that since point[2] can be negative then 
  // the rounding off is sometimes 'greater' than point[2]
  G4int iz = int(z + cmsExp::noffZ);

  // The previous implementation was:
  //    G4double dz = point[2]- int(point[2]);
  // the main issue is that sicne point[2] can be negative then,
  // dz can be negative (which is fatal for the interpolation!
  G4double dz = z + cmsExp::noffZ - iz;
  
  G4double Bz_lb = fieldMap[iz][ir].x();
  G4double Bz_ub = fieldMap[iz+1][ir].x();

  bField[2] = tesla*(Bz_lb + (Bz_ub-Bz_lb)*dz);  
  
  if(rho>0) {
    G4double Br_lb = fieldMap[iz][ir].y();
    G4double Br_ub = fieldMap[iz][ir+1].y();
    bField[0] = tesla*(Br_lb + (Br_ub-Br_lb)*dr)*point[0]/rho;
    bField[1] = tesla*(Br_lb + (Br_ub-Br_lb)*dr)*point[1]/rho;
  }
}

void cmsExpMagneticField::GetParametricBfield(G4double const *point, 
			              G4double *bField, G4double rho)  const 
{
  bool isTracker = (std::abs(point[2])< cmsExp::trackerZmax && 
		    rho < cmsExp::trackerRmax ) ? true: false ;

  if (isTracker) {
    //    theTrakcerBfield->getBxyz(point, bField);
  }
  else {
    //use the volumebase magnetic field outside the tracker
    GetVolumeBaseBfield(point, bField);
  }
}

void cmsExpMagneticField::GetUniformBfield(G4double const *point, 
			              G4double *bField, G4double rho)  const 
{
  bool isTracker = (std::abs(point[2])< cmsExp::trackerZmax && 
		    rho < cmsExp::trackerRmax ) ? true: false ;

  if (isTracker) {
    bField[2] = flux;
  }
  else {
    //temporarily use the volumebase magnetic field outside the tracker 
    GetVolumeBaseBfield(point, bField);
  }

}
