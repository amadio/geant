#include <TFile.h>
#include <TGeoExtension.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TList.h>
#include <TRandom.h>
#include <TGeoBBox.h>
#include <GeantTrack.h>

#include <TPartIndex.h>
#include <TMXsec.h>
#include <TPXsec.h>

#include "base/Global.h"
using vecgeom::kTwoPi;

void GenerateEvent(double avemult, double energy, double fVertex[3]);
double SampleMaxwell(double emean);
void IncreaseStack();
void VertexIn(TGeoBBox *bbox, double ori[3]);

using namespace Geant;

static int stacksize = 100;
static int hwmark = 0;
static GeantTrack *particleStack = new GeantTrack[stacksize];
static TGeoManager *geom;

int main(int argc, char *argv[]) {

  for (int i = 0; i < argc; ++i) {
    printf("argv[%d] = %s\n", i, argv[i]);
  }
  /*
     int nevent=1;
     if(argc>1) sscanf(argv[1],"%d",&nevent);

     double avemult = 10.;
     if(argc>2) sscanf(argv[2],"%lf",&avemult);

     double energy = 10.;
     if(argc>3) sscanf(argv[3],"%lf",&energy);

     printf("Generating %d events with ave multiplicity %f and energy %f\n",nevent,avemult, energy);

     const char *geofile="http://root.cern.ch/files/cms.root";
     geom = TGeoManager::Import(geofile);

     // loop materials
     TFile *f = new TFile("xsec.root");
     f->Get("PartIndex");

     TPXsec::SetVerbose(1);
     TPartIndex::I()->SetEnergyGrid(1e-3,1e3,100);
     TList *matlist = (TList*) geom->GetListOfMaterials();
     TIter next(matlist);
     TGeoMaterial *mat=0;

     TList *matXS = new TList();
     matXS->Add(TPartIndex::I());
     TMXsec *mxs=0;
     printf("Total of %d materials\n",matlist->GetSize());
     while((mat = (TGeoMaterial*) next())) {
        if(!mat->IsUsed()) continue;
        int nelem = mat->GetNelements();
        int *z = new int[nelem];
        int *a = new int[nelem];
        float *w = new float[nelem];
        for(int iel=0; iel<nelem; ++iel) {
           double ad;
           double zd;
           double wd;
           mat->GetElementProp(ad,zd,wd,iel);
           a[iel]=ad;
           z[iel]=zd;
           w[iel]=wd;
           //	 printf("Mixture %s element %s z %d a %d\n",
           //	mat->GetName(), mat->GetElement(iel)->GetName(),
           //	z[iel],a[iel]);
        }
        mxs = new TMXsec(mat->GetName(),mat->GetTitle(),
                         z,a,w,nelem,mat->GetDensity(),kTRUE);
        matXS->Add(mxs);
        mat->SetFWExtension(new TGeoRCExtension(mxs));
        //      myObject = mat->GetExtension()->GetUserObject();
        delete [] a;
        delete [] z;
        delete [] w;
     }

     TMXsec::Prune();

     TFile *fmxs = new TFile("mxs.root","recreate");
     fmxs->SetCompressionLevel(0);
     matXS->Write();
     fmxs->Close();

     TGeoVolume *top = geom->GetTopVolume();
     TGeoShape *shape = top->GetShape();
     TGeoBBox *bbox = (TGeoBBox*) shape;
     double dx = bbox->GetDX();
     double dy = bbox->GetDY();
     double dz = bbox->GetDZ();
     const double *origin = bbox->GetOrigin();
     printf("Top volume is %s shape %s\n",top->GetName(),shape->GetName());
     printf("BBox dx %f dy %f dz %f origin %f %f %f\n",dx,dy,dz,origin[0],origin[1],origin[2]);

     for(int iev=0; iev<nevent; ++iev) {
        // should define a vertex, origin for the moment
        double vertex[3]={0,0,0};
        VertexIn(bbox,vertex);
        GenerateEvent(avemult, energy, vertex);
        while(hwmark) {
           GeantTrack *track = &particleStack[--hwmark];
           printf("Transporting particle #%d %s energy %g\n",hwmark+1,
                  TDatabasePDG::Instance()->GetParticle(track->fPDG)->GetName(),
                  track->fE-track->fMass);
           double pos[3] = {track->fXpos,track->fYpos,track->fZpos};
           double dir[3] = {track->fXdir,track->fYdir,track->fZdir};
           TGeoNode *current = geom->InitTrack(pos,dir);
           double pintl = -log(gRandom->Rndm());
           printf("Initial point %f %f %f in %s\n",pos[0],pos[1],pos[2],current->GetName());
           while(!geom->IsOutside()) {
              mat = current->GetVolume()->GetMaterial();
              const double *cpos = geom->GetCurrentPoint();
              printf("Point now %f %f %f in %s made of %s\n",cpos[0],cpos[1],cpos[2],
                     current->GetName(),current->GetVolume()->GetMaterial()->GetName());
              double ken = track->fE-track->fMass;
              TMXsec *mx = ((TMXsec *)
                            ((TGeoRCExtension*)
                             mat->GetFWExtension())->GetUserObject());
              double xlen = mx->Xlength(track->fGVcode,ken,sqrt(
     (track->fE+track->fMass)*(track->fE-track->fMass) )  );
              double pnext = pintl*xlen;
              current = geom->FindNextBoundaryAndStep(pnext);
              double snext = geom->GetStep();
              printf("pnext = %f snext = %f\n",pnext,snext);
              if(pnext<=snext) {
                 //phys wins
                 int reac;
                 mat->Print();
                 TEXsec *el = mx->SampleInt(track->fGVcode,ken,reac);
                 printf("particle does a %s on %s\n",TPartIndex::I()->ProcName(reac),el->GetName());
                 break;
              } else {
                 // geom wins
                 pintl-=snext/xlen;
                 //	       printf("Geom wins\n");
              }
              //	    dirnew = something;
              //	       geom->SetCurrentDirection(dirnew);
           }
           if(geom->IsOutside()) printf("Particle exited setup\n");
        }
     }
  */
  /*
  double dir[3];
  double pos[3];
  TGeoNode *current=0;
  TGeoNode *nexnode=0;
  TIter next(particleStack);
  for(int iev=0; iev<nevent; ++iev) {
     GenerateEvent();
     next.Reset();
     GeantTrack *tr=0;
     while((tr=(GeantTrack*)next())) {
        int GVindex = TPartIndex::I()->PartIndex(tr->pdg);
        tr->Direction(dir);
        x[0]=tr->xpos;
        x[1]=tr->ypos;
        x[2]=tr->zpos;
        // where am I

     }
  */
  return 0;
}

#define NPART 11

void GenerateEvent(double avemult, double energy, double fVertex[3]) {
  static bool first = kTRUE;
  static const int kMaxPart = NPART;
  static const char *GVname[NPART] = {"pi+", "pi-", "proton", "antiproton", "neutron", "antineutron",
                                      "e-",  "e+",  "gamma",  "mu+",        "mu-"};
  static const Species_t GVspecies[NPART] = {kHadron, kHadron, kHadron, kHadron, kHadron, kHadron,
                                             kLepton, kLepton, kLepton, kLepton, kLepton};
  static int GVpart[NPART];
  static float GVprob[NPART] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};

  const double etamin = -3, etamax = 3;

  // Initialise simple generator
  if (first) {
    double sumprob = 0;
    for (int ip = 0; ip < kMaxPart; ++ip) {
      GVpart[ip] = TPartIndex::I()->PartIndex(GVname[ip]);
      printf("part %s code %d\n", GVname[ip], GVpart[ip]);
      sumprob += GVprob[ip];
    }
    for (int ip = 0; ip < kMaxPart; ++ip) {
      GVprob[ip] /= sumprob;
      if (ip)
        GVprob[ip] += GVprob[ip - 1];
    }
    first = kFALSE;
  }

  int ntracks = gRandom->Poisson(avemult) + 0.5;

  hwmark = 0;
  for (int i = 0; i < ntracks; i++) {
    if (hwmark == stacksize)
      IncreaseStack();
    GeantTrack *track = &particleStack[hwmark++];
    double prob = gRandom->Uniform();
    for (int j = 0; j < kMaxPart; ++j) {
      if (prob <= GVprob[j]) {
        track->fGVcode = GVpart[j];
        track->fPDG = TPartIndex::I()->PDG(GVpart[j]);
        track->fSpecies = GVspecies[j];
#ifdef USE_VECGEOM_NAVIGATOR
        printf("Generating a %s\n", Particle::GetParticle(track->fPDG).Name());
#else
        printf("Generating a %s\n", TDatabasePDG::Instance()->GetParticle(track->fPDG)->GetName());
#endif
        //	    pdgCount[j]++;
        break;
      }
    }
    if (!track->fPDG)
#ifdef USE_VECGEOM_NAVIGATOR
    {
      std::cout << __func__ << ":: No Particle generated!" << std::endl;
      exit(1);
    }
    const Particle *part = &Particle::GetParticle(track->fPDG);
    track->fCharge = part->Charge();
#else
      Fatal("ImportTracks", "No particle generated!");
    TParticlePDG *part = TDatabasePDG::Instance()->GetParticle(track->fPDG);
    track->fCharge = part->Charge() / 3.;
#endif
    track->fMass = part->Mass();
    track->fXpos = fVertex[0];
    track->fYpos = fVertex[1];
    track->fZpos = fVertex[2];
    double ekin = SampleMaxwell(energy / avemult);
    track->fE = ekin + track->fMass;
    double p = sqrt(ekin * (2 * ekin + track->fMass));
    double eta = gRandom->Uniform(etamin, etamax); // multiplicity is flat in rapidity
    double theta = 2 * atan(exp(-eta));
    // double theta = acos((1.-2.*gRandom->Rndm()));
    double phi = kTwoPi * gRandom->Rndm();
    track->fP = p;
    track->fXdir = sin(theta) * cos(phi);
    track->fYdir = sin(theta) * sin(phi);
    track->fZdir = cos(theta);
    track->fFrombdr = kFALSE;
  }
  //      Printf("Event #%d: Generated species for %6d particles:", event, ntracks);
}

double SampleMaxwell(double emean) {
  double th = gRandom->Uniform() * kTwoPi;
  double rho = sqrt(-log(gRandom->Uniform()));
  double mx = rho * sin(th);
  return 2 * emean * (-log(gRandom->Uniform()) + mx * mx) / 3.;
}

void IncreaseStack() {
  int newstacksize = stacksize * 1.5;
  GeantTrack *tmp = new GeantTrack[newstacksize];
  memcpy((void *)tmp, (void *)particleStack, stacksize * sizeof(GeantTrack));
  delete[] particleStack;
  particleStack = tmp;
  stacksize = newstacksize;
}

void VertexIn(TGeoBBox *bbox, double ori[3]) {
  double eta = 0;
  do {
    eta = gRandom->Rndm();
    ori[0] = bbox->GetDX() * eta * eta * (gRandom->Rndm() > 0.5 ? 1 : -1);
    eta = gRandom->Rndm();
    ori[1] = bbox->GetDY() * eta * eta * (gRandom->Rndm() > 0.5 ? 1 : -1);
    eta = gRandom->Rndm();
    ori[2] = bbox->GetDZ() * eta * eta * (gRandom->Rndm() > 0.5 ? 1 : -1);
    geom->SetCurrentPoint(ori);
  } while (geom->IsOutside());
}
