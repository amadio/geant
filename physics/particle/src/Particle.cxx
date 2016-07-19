#include "Particle.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstring>
#ifndef GEANT_NVCC
using std::cout;
#endif
using std::endl;
using std::ofstream;
using std::ifstream;
using std::stringstream;
using std::setw;
using std::setfill;

static const double kPlankBar = 6.5821192815e-25; // GeV s

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

#ifdef GEANT_NVCC
GEANT_CUDA_DEVICE_CODE vecgeom::map<int, Particle> *fParticlesDev = nullptr;
vecgeom::map<int, Particle> *fParticlesHost                       = nullptr;

GEANT_CUDA_BOTH_CODE
char *strncpy(char *dest, const char *src, size_t n)
{
  char *ret = dest;
  do {
    if (!n--) return ret;
  } while (*dest++ = *src++);
  while (n--)
    *dest++ = 0;
  return ret;
}

#else
std::map<int, Particle> *Particle::fParticles=0;
#endif

#ifndef GEANT_NVCC
std::ostream &operator<<(std::ostream &os, const Particle &part)
{
  os << part.fName << "(" << part.fPDG << ") Class:" << part.fClass << " Q:" << part.fCharge << " m:" << part.fMass
     << " lt:" << part.fLife << " I:" << (int)part.fIsospin << " I3:" << (int)part.fIso3 << " S:" << (int)part.fStrange
     << " F:" << (int)part.fFlavor << " #:" << (int)part.fNdecay << " code:" << (int)part.fCode << endl;
  return os;
}
#endif
//________________________________________________________________________________________________

GEANT_CUDA_BOTH_CODE
Particle::Particle()
    : fPDG(0), fMatter(true), fPcode(0), fCharge(0), fMass(-1), fWidth(0), fIsospin(0), fIso3(0), fStrange(0),
      fFlavor(0), fTrack(0), fNdecay(0), fCode(-1)
{
  strncpy(fName, "Default", 8);
  strncpy(fClass, "", 1);
}

//________________________________________________________________________________________________
GEANT_CUDA_BOTH_CODE
Particle::Particle(const char *name, int pdg, bool matter, const char *pclass, int pcode, double charge, double mass,
                   double width, int isospin, int iso3, int strange, int flavor, int track, int code)
    : fPDG(pdg), fMatter(matter), fPcode(pcode), fCharge(charge), fMass(mass), fWidth(width), fIsospin(isospin),
      fIso3(iso3), fStrange(strange), fFlavor(flavor), fTrack(track), fNdecay(0), fCode(code)
{
  strncpy(fName, name, 255);
  fName[255] = '\0';
  strncpy(fClass, pclass, 255);
  fClass[255] = '\0';

#ifndef GEANT_NVCC
  if (!fParticles) fParticles = new Map_t;
  if (fParticles->count(fPDG) != 0) {
    printf("Particle %d already there\n", fPDG);
    return;
  }

  if (fWidth > 0)
    fLife = kPlankBar / fWidth;
  else
    fLife             = 0;
  (*fParticles)[fPDG] = *this;
#else
#ifndef GEANT_CUDA_DEVICE_CODE
  if (!fParticlesHost) fParticlesHost = new Map_t;
  if (fParticlesHost->count(fPDG) != 0) {
    printf("Particle %d already there\n", fPDG);
    return;
  }

  if (fWidth > 0)
    fLife = kPlankBar / fWidth;
  else
    fLife                 = 0;
  (*fParticlesHost)[fPDG] = *this;

#else
  if (!fParticlesDev) fParticlesDev = new Map_t;
  if (fParticlesDev->count(fPDG) != 0) {
    printf("Particle %d already there\n", fPDG);
    return;
  }

  if (fWidth > 0)
    fLife = kPlankBar / fWidth;
  else
    fLife                = 0;
  (*fParticlesDev)[fPDG] = *this;
#endif
#endif
}
GEANT_CUDA_BOTH_CODE
Particle::Particle(const Particle &other)
    : fPDG(other.fPDG), fMatter(other.fMatter), fPcode(other.fPcode), fCharge(other.fCharge), fMass(other.fMass),
      fWidth(other.fWidth), fIsospin(other.fIsospin), fStrange(other.fStrange), fFlavor(other.fFlavor),
      fTrack(other.fTrack), fCode(other.fCode)
{

  strncpy(fName, other.fName, 255);
  fName[255] = '\0';
  strncpy(fClass, other.fClass, 255);
  fClass[255] = '\0';
}
GEANT_CUDA_BOTH_CODE
Particle &Particle::operator=(const Particle &part)
{
  if (&part != this) {
    strncpy(fName, part.fName, 255);
    fName[255] = '\0';
    strncpy(fClass, part.fClass, 255);
    fClass[255] = '\0';
    fPDG        = part.fPDG;
    fMatter     = part.fMatter;
    fPcode      = part.fPcode;
    fCharge     = part.fCharge;
    fMass       = part.fMass;
    fWidth      = part.fWidth;
    fLife       = part.fLife;
    fIsospin    = part.fIsospin;
    fIso3       = part.fIso3;
    fStrange    = part.fStrange;
    fFlavor     = part.fFlavor;
    fTrack      = part.fTrack;
    fCode       = part.fCode;
    for (unsigned int i = 0; i < part.fDecayList.size(); ++i)
      fDecayList[i]     = part.fDecayList[i];
    fNdecay = (unsigned char) fDecayList.size();
    // cout << "Particle::= called" << endl;
  }
  return *this;
}

//________________________________________________________________________________________________
#ifndef GEANT_NVCC
void Particle::ReadFile(std::string infilename, bool output)
{
  int count;
  std::string name;
  int pdg;
  bool matter;
  int pcode;
  std::string pclass;
  int charge;
  double mass, width;
  int isospin, iso3, strange, flavor, track, ndecay;
  int ipart, acode;
  int kcount       = 0;
  const int ksplit = 40;
  int kfunc        = 0;

  std::ofstream outfile;
  ifstream infile(infilename);
  std::string line;
  std::stringstream outline;
  std::stringstream filename;

  int np = 0;
  while (getline(infile, line)) {
    if (np == 0)
      for (int i = 0; i < 3; ++i) {
        if (line.substr(0, 1) != "#") {
          printf("There should be three lines of comments at the beginning \n");
          exit(1);
        }
        getline(infile, line);
      }
    if (line.substr(0, 1) == "#") {
      printf("Disaster embedded comment!!!\n");
      exit(1);
    }

    ++np;
    //      cout << "line: " << line << endl;
    GetPart(line, count, name, pdg, matter, pcode, pclass, charge, mass, width, isospin, iso3, strange, flavor, track,
            ndecay, ipart, acode);
    if (np != count) {
      printf("Disaster count np(%d) != count(%d)", np, count);
      exit(1);
    }
    /*    
    cout << " count:" << count << " name:" << name << " pdg:" << pdg << " matter:" << matter << " pcode:" << pcode
   <<" pclass:" << pclass << " charge:" << charge << " mass:" << mass << " width:" << width << " isospin:"
   << isospin << " iso3:" << iso3 << " strange:" << strange << " flavor:" << flavor << " track:" << track
   << " ndecay:" << ndecay << " ipart:" << ipart << " acode:" << acode << endl;
    */
    if (pdg >= 0) {
      if (isospin != -100) isospin = (isospin - 1) / 2;
      if (strange != -100) strange = (strange - 1) / 2;
      new Particle(name.c_str(), pdg, matter, pclass.c_str(), pcode, charge / 3., mass, width, isospin, iso3, strange,
                   flavor, track);
      Particle &part = fParticles->at(pdg);
      if (ndecay > 0) {
        for (int i = 0; i < 3; ++i) {
          getline(infile, line);
          if (line.substr(0, 1) != "#") {
            printf("Disaster comment!!!\n");
            exit(1);
          }
        }
        VectorDecay_t decaylist;
        decaylist.clear();
        for (int i = 0; i < ndecay; ++i) {
          getline(infile, line);
          if (line.substr(0, 1) == "#") {
            printf("Disaster embedded comment!!!\n");
            exit(1);
          }
          Decay decay;
          int dcount;
          GetDecay(line, dcount, decay);
          if (dcount != i + 1) {
            printf("dcount (%d) != i+1 (%d)", dcount, i + 1);
            exit(1);
          }
	  // cout << "         " << dcount << " " << decay << endl;
          part.AddDecay(decay);
        }
        part.NormDecay();
      }
    } else {
      // Create antiparticle from particle
      if (fParticles->find(-pdg) == fParticles->end()) {
        printf("Cannot build the antiparticle because the particle %d is not there!", -pdg);
        exit(1);
      }
      Particle &p = fParticles->at(-pdg);
      new Particle(name.c_str(), pdg, matter, p.Class(), p.Pcode(), p.Charge() == 0 ? 0 : -p.Charge(), p.Mass(),
                   p.Width(), p.Isospin() == 0 ? 0 : -p.Isospin(), p.Isospin() == 0 ? 0 : -p.Isospin(),
                   p.Iso3() == 0 ? 0 : -p.Iso3(), p.Strange() == 0 ? 0 : -p.Strange(),
                   p.Flavor() == 0 ? 0 : -p.Flavor(), p.Track());
      Particle &part = fParticles->at(pdg);
      Decay decay;
      auto dl = p.DecayList();
      for (int i = 0; i < p.Ndecay(); ++i) {
        decay.Clear();
        decay.SetType(dl.at(i).Type());
        decay.SetBr(dl.at(i).Br());
        for (int j = 0; j < dl.at(i).NDaughters(); ++j)
          decay.AddDaughter(-dl.at(i).Daughter(j));
        //	    cout << "         " << i << " " << decay << endl;
        part.AddDecay(decay);
      }
    }
  }

  if (output) {

    bool partdef = false;
    for (auto p = Particle::Particles().begin(); p != Particle::Particles().end(); ++p) {
      if (kcount % ksplit == 0) {
	 if (kcount > 0) {
	    outline << "}" << endl << endl;
	    outline << " } // End of inline namespace" << endl;
	    outline << " } // End of geant namespace" << endl;
	    outline << "#if defined(__clang__) && !defined(__APPLE__)" << endl;
	    outline << "#pragma clang optimize on" << endl;
	    outline << "#endif" << endl;
	    filename << "CreateParticle" << setfill('0') << setw(4) << kfunc-1 << ".cxx" 
		     << setfill(' ') << setw(0);
	    outfile.open(filename.str());
	    outfile << outline.str();
	    outfile.close();
	    filename.str("");
	    filename.clear();
	    outline.str("");
	    outline.clear(); // Clear state flags.
	 }

	    
	 outline << "#if defined(__clang__) && !defined(__APPLE__)" << endl;
	 outline << "#pragma clang optimize off" << endl;
	 outline << "#endif" << endl;
	 outline << "#include \"Particle.h\"" << endl;
	 outline << "#ifdef GEANT_NVCC" << endl << "#include \"base/Vector.h\"" << endl
		 << "template <typename T>" << endl << "using vector = vecgeom::Vector<T>;" << endl
		 << "#else" << endl << "using std::vector;" << endl << "#endif" << endl;
	 outline << "namespace geant {" << endl;
	 outline << "   inline namespace GEANT_IMPL_NAMESPACE {" << endl << endl;
	 outline << endl << "//" << setw(80) << setfill('_') << "_" << endl << setfill(' ') << setw(0);
	 outline << "GEANT_CUDA_BOTH_CODE" << endl << "void CreateParticle" << setfill('0') << setw(4) << kfunc << "() {" << endl
                << setfill(' ') << setw(0);
        partdef = false;
        ++kfunc;
      }
      ++kcount;
      const Particle &part = p->second;
      std::string name(part.Name());
      outline << endl << "   // Creating " << name << endl;
      size_t quote = name.find_first_of("\\\"");
      if (quote != std::string::npos)
        name = name.substr(0, quote) + "\\\"" + name.substr(quote + 1, name.size() - quote);
      outline << "   new Particle("
              << "\"" << name << "\", " << part.PDG() << ", " << part.Matter() << ", "
              << "\"" << part.Class() << "\", " << part.Pcode() << ", " << part.Charge() << ", " << part.Mass() << ", "
              << part.Width() << ", " << part.Isospin() << ", " << part.Iso3() << ", " << part.Strange() << ", "
              << part.Flavor() << ", " << part.Track() << ");" << endl;
      //	 cout << name << " " << part.Ndecay() << endl;
      if (part.Ndecay() > 0) {
        if (!partdef) {
          outline << "   Particle *part = 0;" << endl;
          partdef = true;
        }
        outline << "   part = const_cast<Particle*>(&Particle::Particles().at(" << part.PDG() << "));" << endl;
        for (VectorDecay_t::const_iterator d = part.DecayList().begin(); d != part.DecayList().end(); ++d) {
          outline << "   part->AddDecay(Particle::Decay(" << d->Type() << ", " << d->Br() << ", "
                  << " vector<int>{";
          for (int i = 0; i < d->NDaughters(); ++i) {
            outline << d->Daughter(i);
            if (i < d->NDaughters() - 1) outline << ",";
          }
          outline << "}));" << endl;
        }
      }
    }

    if (kcount % ksplit != 0) {
       outline << "}" << endl << endl;
       outline << " } // End of inline namespace" << endl;
       outline << " } // End of geant namespace" << endl;
       outline << "#if defined(__clang__) && !defined(__APPLE__)" << endl;
       outline << "#pragma clang optimize on" << endl;
       outline << "#endif" << endl;
       filename << "CreateParticle" << setfill('0') << setw(4) << kfunc-1 << ".cxx" 
			     << setfill(' ') << setw(0);
       outfile.open(filename.str());
       outfile << outline.str();
       outfile.close();
       filename.str("");
       filename.clear();
       outline.str("");
       outline.clear(); // Clear state flags.
    }

    outfile.open("CreateParticles.cxx");
    outfile << "#include \"Particle.h\"" << endl;
    outfile << "namespace geant {" << endl;
    outfile << "   inline namespace GEANT_IMPL_NAMESPACE {" << endl << endl;
    for (int i = 0; i < kfunc; ++i)
      outfile << "void CreateParticle" << setfill('0') << setw(4) << i << "();" << endl;

    
    outfile << endl << "//" << setw(80) << setfill('_') << "_" << endl << setfill(' ') << setw(0);
    outfile << "GEANT_CUDA_BOTH_CODE" << endl << "void Particle::CreateParticles() {" << endl;
    outfile << "   static bool initDone=false;" << endl;
    outfile << "   if(initDone) return;" << endl;
    outfile << "   initDone = true;" << endl;
    for (int i = 0; i < kfunc; ++i)
      outfile << "    CreateParticle" << setfill('0') << setw(4) << i << "();" << endl;
    outfile << "}" << endl;
    outfile << " } // End of inline namespace" << endl;
    outfile << " } // End of geant namespace" << endl;
    outfile << "#if defined(__clang__) && !defined(__APPLE__)" << endl;
    outfile << "#pragma clang optimize on" << endl;
    outfile << "#endif" << endl;
    outfile.close();
  }
}
//________________________________________________________________________________________________
void Particle::GetDecay(const std::string &line, int &dcount, Decay &decay)
{
  int dtype;
  double br;
  int ndec;
  int daughter;

  std::stringstream ss(line);
  decay.Clear();
  ss >> dcount;
  ss >> dtype;
  ss >> br;
  decay.SetType(dtype);
  decay.SetBr(br);
  ss >> ndec;
  for (int i = 0; i < ndec; ++i) {
    ss >> daughter;
    decay.AddDaughter(daughter);
  }
}

//________________________________________________________________________________________________
void Particle::NormDecay()
{
  double brt = 0;

  int ndec = fDecayList.size();
  if (ndec) {
    for (VectorDecay_t::iterator idec = fDecayList.begin(); idec != fDecayList.end(); ++idec)
      brt += idec->Br();
    if (brt) {
      brt = 1 / brt;
      for (VectorDecay_t::iterator idec = fDecayList.begin(); idec != fDecayList.end(); ++idec)
        idec->SetBr(idec->Br() * brt);
      for (unsigned int i = 0; i < fDecayList.size() - 1; ++i)
        for (unsigned int j = i + 1; j < fDecayList.size(); ++j)
          if (fDecayList.at(i).Br() < fDecayList.at(j).Br()) {
            Decay dec        = fDecayList.at(i);
            fDecayList.at(i) = fDecayList.at(j);
            fDecayList.at(j) = dec;
          }
    }
  } else if (fLife == 0)
    fLife = 1e38;
}

//________________________________________________________________________________________________
void Particle::GetPart(const std::string &line, int &count, std::string &name, int &pdg, bool &matter, int &pcode,
                       std::string &pclass, int &charge, double &mass, double &width, int &isospin, int &iso3,
                       int &strange, int &flavor, int &track, int &ndecay, int &ipart, int &acode)
{

  count   = 0;
  name    = "";
  pdg     = 0;
  matter  = false;
  pclass  = "";
  charge  = 0;
  mass    = 0;
  width   = 0;
  isospin = 0;
  iso3    = 0;
  strange = 0;
  flavor  = 0;
  track   = 0;
  ndecay  = 0;
  ipart   = 0;
  acode   = 0;
  pcode   = 0;

  stringstream ss(line);

  ss >> count;
  ss >> name;
  ss >> pdg;

  if (pdg < 0) {
    matter = false;
    ss >> ipart;
    ss >> acode;
  } else {

    int imatter = 0;
    ss >> imatter;
    matter = (imatter == 1);

    ss >> pcode;
    ss >> pclass;
    ss >> charge;
    ss >> mass;
    ss >> width;
    ss >> isospin;
    ss >> iso3;
    ss >> strange;
    ss >> flavor;
    ss >> track;
    ss >> ndecay;
  }
}
std::ostream &operator<<(std::ostream &os, const Particle::Decay &dec)
{
  os << "Type " << static_cast<int>(dec.fType) << " br " << dec.fBr << " products: ";
  for (unsigned int i = 0; i < dec.fDaughters.size(); ++i)
    os << " " << dec.fDaughters[i];

  /*   for(int i=0; i<dec.fDaughters.size(); ++i)
 os << " " << Particle::Particles().at(dec.fDaughters[i]).Name(); */

  return os;
}
#endif
}
}
