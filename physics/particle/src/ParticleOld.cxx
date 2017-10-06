#include "ParticleOld.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstring>
#ifndef VECCORE_CUDA
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

#ifdef VECCORE_CUDA
VECCORE_ATT_DEVICE vecgeom::map<int, ParticleOld> *fParticlesDev = nullptr;
vecgeom::map<int, ParticleOld> *fParticlesHost                       = nullptr;

VECCORE_ATT_HOST_DEVICE
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
std::map<int, ParticleOld> *ParticleOld::fParticles=0;
#endif

#ifndef VECCORE_CUDA
std::ostream &operator<<(std::ostream &os, const ParticleOld &part)
{
  os << part.fName << "(" << part.fPDG << ") Class:" << part.fClass << " Q:" << part.fCharge << " m:" << part.fMass
     << " lt:" << part.fLife << " I:" << (int)part.fIsospin << " I3:" << (int)part.fIso3 << " S:" << (int)part.fStrange
     << " F:" << (int)part.fFlavor << " #:" << (int)part.fNdecay << " code:" << (int)part.fCode << endl;
  return os;
}
#endif

//________________________________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
ParticleOld::ParticleOld()
    : fPDG(0), fMatter(true), fPcode(0), fCharge(0), fMass(-1), fWidth(0), fIsospin(0), fIso3(0), fStrange(0),
      fFlavor(0), fTrack(0), fNdecay(0), fCode(-1)
{
  strncpy(fName, "Default", 8);
  strncpy(fClass, "", 1);
}

//________________________________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
ParticleOld::ParticleOld(const char *name, int pdg, bool matter, const char *pclass, int pcode, double charge, double mass,
                   double width, int isospin, int iso3, int strange, int flavor, int track, int code)
    : fPDG(pdg), fMatter(matter), fPcode(pcode), fCharge(charge), fMass(mass), fWidth(width), fIsospin(isospin),
      fIso3(iso3), fStrange(strange), fFlavor(flavor), fTrack(track), fNdecay(0), fCode(code)
{
  strncpy(fName, name, 255);
  fName[255] = '\0';
  strncpy(fClass, pclass, 255);
  fClass[255] = '\0';

#ifndef VECCORE_CUDA
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
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
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

VECCORE_ATT_HOST_DEVICE
ParticleOld::ParticleOld(const ParticleOld &other)
    : fPDG(other.fPDG), fMatter(other.fMatter), fPcode(other.fPcode), fCharge(other.fCharge), fMass(other.fMass),
      fWidth(other.fWidth), fIsospin(other.fIsospin), fStrange(other.fStrange), fFlavor(other.fFlavor),
      fTrack(other.fTrack), fCode(other.fCode)
{

  strncpy(fName, other.fName, 255);
  fName[255] = '\0';
  strncpy(fClass, other.fClass, 255);
  fClass[255] = '\0';
}
VECCORE_ATT_HOST_DEVICE
ParticleOld &ParticleOld::operator=(const ParticleOld &part)
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
    // cout << "ParticleOld::= called" << endl;
  }
  return *this;
}

#ifndef VECCORE_CUDA

// Returns the number of lines written.
static unsigned WriteOneParticle(std::stringstream &outline, const ParticleOld &part, unsigned int nfuncs)
{
  unsigned int nlines = 0;
  std::string name(part.Name());
  if (nfuncs) {
    outline << "VECCORE_ATT_HOST_DEVICE" << endl << "static void InternalCreateParticle" << setfill('0')
            << setw(4) << nfuncs-1 << "() {" << endl
            << setfill(' ') << setw(0);
    if (part.Ndecay() > 0) {
      outline << "   vector<int> daughters;\n";
      outline << "   ParticleOld *part = nullptr;\n";
    }
  }
  outline << endl << "   // Creating " << name << endl;
  size_t quote = name.find_first_of("\\\"");
  if (quote != std::string::npos)
    name = name.substr(0, quote) + "\\\"" + name.substr(quote + 1, name.size() - quote);
  outline << "   new ParticleOld("
          << "\"" << name << "\", " << part.PDG() << ", " << part.Matter() << ", "
          << "\"" << part.Class() << "\", " << part.Pcode() << ", " << part.Charge() << ", " << part.Mass() << ", "
          << part.Width() << ", " << part.Isospin() << ", " << part.Iso3() << ", " << part.Strange() << ", "
          << part.Flavor() << ", " << part.Track() << ");" << endl;
  //	 cout << name << " " << part.Ndecay() << endl;
  ++nlines;
  if (part.Ndecay() > 0) {
    outline << "   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(" << part.PDG() << "));" << endl;
    for (ParticleOld::VectorDecay_t::const_iterator d = part.DecayList().begin(); d != part.DecayList().end(); ++d) {
      outline << "   daughters.clear();\n";
      ++nlines;
      for (int i = 0; i < d->NDaughters(); ++i) {
        outline << "   daughters.push_back(" << d->Daughter(i) << ");\n";
        ++nlines;
      }
      outline << "   part->AddDecay(ParticleOld::Decay(" << d->Type() << ", " << d->Br() << ", "
              << " daughters));\n";
      ++nlines;
    }
  }
  if (nfuncs) {
    outline << "}\n\n";
  }
  return nlines;
}

static void SwitchFile(std::stringstream &outline, unsigned int nfiles, unsigned int nfunc = 0)
{
  if (nfunc) {
    outline << "} // anonymous namespace\n\n";
    outline << "VECCORE_ATT_HOST_DEVICE" << endl << "void CreateParticle" << setfill('0') << setw(4) << nfiles-1 << "() {" << endl
            << setfill(' ') << setw(0);
    for (unsigned int i = 0; i < nfunc; ++i)
      outline << "  InternalCreateParticle" << setfill('0') << setw(4) << i << "();" << endl;
  }
  outline << "}" << endl << endl;
  outline << "} // End of inline namespace" << endl;
  outline << "} // End of geant namespace" << endl;
  outline << "#if defined(__clang__) && !defined(__APPLE__)" << endl;
  outline << "#pragma clang optimize on" << endl;
  outline << "#endif" << endl;

  std::stringstream filename;
  filename << "CreateParticle" << setfill('0') << setw(4) << nfiles-1 << ".cxx"
           << setfill(' ') << setw(0);

  std::ofstream outfile;
  outfile.open(filename.str());
  outfile << outline.str();
  outfile.close();

  outline.str("");
  outline.clear(); // Clear state flags.
}

static void StartFile(std::stringstream &outline, unsigned int nfiles, bool oneFuncPer = false)
{
  outline << "// This files was autogenerated by geant::ParticleOld::ReadFile\n\n";
  outline << "#if defined(__clang__) && !defined(__APPLE__)" << endl;
  outline << "#pragma clang optimize off" << endl;
  outline << "#endif" << endl;
  outline << "#include \"ParticleOld.h\"" << endl;
  outline << "#ifdef VECCORE_CUDA" << endl << "#include \"base/Vector.h\"" << endl
          << "template <typename T>" << endl << "using vector = vecgeom::Vector<T>;" << endl
          << "#else" << endl << "using std::vector;" << endl << "#endif" << endl << endl;
  outline << "namespace geant {" << endl;
  outline << "inline namespace GEANT_IMPL_NAMESPACE {" << endl << endl;
  outline << endl << "//" << setw(80) << setfill('_') << "_" << endl << setfill(' ') << setw(0);
  if (!oneFuncPer) {
    outline << "VECCORE_ATT_HOST_DEVICE" << endl << "void CreateParticle" << setfill('0') << setw(4) << nfiles << "() {" << endl
            << setfill(' ') << setw(0);
    outline << "   vector<int> daughters;\n";
    outline << "   ParticleOld *part = nullptr;\n";
  } else {
    outline << "namespace {\n";
  }
}

//________________________________________________________________________________________________
void ParticleOld::ReadFile(std::string infilename, bool output)
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
  int nparticles = 0;
  // const int ksplit = 20;
  constexpr bool kOneFuncPer = false;
  constexpr unsigned int kMaxLines = 250;
  // On MacOS, 2.7 GHz Intel Core i7
  // To compile and link all the CreateParticle with nvcc
  // number of lines per files
  // 100            -> real   3m58.662s   user   23m12.182s
  // 50 leads to crash in cicc in file 0064.
  // 40 push_back   -> real   2m27.712s   user   16m3.078s
  // 50 pb          -> real   7m6.049s    user   43m44.287s
  // 40 pb          -> real   7m32.432s   user   47m10.548s
  // 100 pb         -> real   5m8.113s    user   32m6.331s
  // 200 pb         -> real   4m0.983s user   26m8.641s
  // 200 pb, oneper -> cicc crash in 0002
  // 225 pb         -> real   3m55.000s   user   25m12.065s    real 4m3.500s  user 27m2.379s
  // 250 pb         -> real   3m48.694s   user   25m24.342s    real 3m42.044s user 24m18.884s
  // 250 pb, oneper -> cicc crash in 0002
  // 275 pb         -> real   4m11.236s   user   25m39.943s    real 3m56.294s user 27m7.945s
  // 300 pb         -> real   3m57.976s   user   24m22.626s
  // 300 pb, oneper -> real   3m52.933s   user   25m21.526s
  // 400 pb         -> real   4m12.210s   user   26m19.758s
  // 400 pb, oneper -> cicc crash in 0001
  // 500 pb         -> real   4m20.674s   user   27m37.601s
  // 500 pb, oneper -> cicc crash in 0001
  // 600 pb         -> real   5m22.091s   user   30m40.852s
  // 600 pb, oneper -> real   4m25.414s   user   26m31.348s


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
      new ParticleOld(name.c_str(), pdg, matter, pclass.c_str(), pcode, charge / 3., mass, width, isospin, iso3, strange,
                   flavor, track);
      ParticleOld &part = fParticles->at(pdg);
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
      ParticleOld &p = fParticles->at(-pdg);
      new ParticleOld(name.c_str(), pdg, matter, p.Class(), p.Pcode(), p.Charge() == 0 ? 0 : -p.Charge(), p.Mass(),
                   p.Width(), p.Isospin() == 0 ? 0 : -p.Isospin(), p.Isospin() == 0 ? 0 : -p.Isospin(),
                   p.Iso3() == 0 ? 0 : -p.Iso3(), p.Strange() == 0 ? 0 : -p.Strange(),
                   p.Flavor() == 0 ? 0 : -p.Flavor(), p.Track());
      ParticleOld &part = fParticles->at(pdg);
      Decay decay;
      auto dl = p.DecayList();
      for (int i = 0; i < p.Ndecay(); ++i) {
        decay.Clear();
        decay.SetType(dl.at(i).Type());
        decay.SetBr(dl.at(i).Br());
        for (int j = 0; j < dl.at(i).NDaughters(); ++j)
          decay.AddDaughter(-dl.at(i).Daughter(j));
        //       cout << "         " << i << " " << decay << endl;
        part.AddDecay(decay);
      }
    }
  }

  if (output) {

    // Approximate number of lines in the current function.
    unsigned int nlines = 0;
    unsigned int nfiles = 0;
    unsigned int nfuncs = 0;
    for (auto p = ParticleOld::Particles().begin(); p != ParticleOld::Particles().end(); ++p) {
      if (nlines == 0 || nlines > kMaxLines) {
        //fprintf(stderr,"Switching to files number %d after count=%d and nlines=%d\n",nfiles,nparticles,nlines);
        if (nparticles > 0) {
          SwitchFile(outline, nfiles, nfuncs);
        }

        StartFile(outline, nfiles, kOneFuncPer);
        ++nfiles;

        nlines = 0;
        nfuncs = 0;
      }
      ++nparticles;
      if (kOneFuncPer) ++nfuncs;
      nlines += WriteOneParticle(outline, p->second, nfuncs);
    }

    if (nlines != 0) {
      SwitchFile(outline, nfiles, nfuncs);
    }

    outfile.open("CreateParticles.cxx");
    outfile << "// This files was autogenerated by geant::ParticleOld::ReadFile\n\n";
    outfile << "#include \"ParticleOld.h\"" << endl;
    outfile << "namespace geant {" << endl;
    outfile << "inline namespace GEANT_IMPL_NAMESPACE {" << endl << endl;
    for (unsigned int i = 0; i < nfiles; ++i) {
      outfile << "VECCORE_ATT_HOST_DEVICE\n";
      outfile << "void CreateParticle" << setfill('0') << setw(4) << i << "();" << endl;
    }

    outfile << "\n#ifdef VECCORE_CUDA\n";
    outfile << "VECCORE_ATT_DEVICE bool fgCreateParticlesInitDoneDev = false;\n";
    outfile << "#endif\n";

    outfile << endl << "//" << setw(80) << setfill('_') << "_" << endl << setfill(' ') << setw(0);
    outfile << "VECCORE_ATT_HOST_DEVICE" << endl << "void ParticleOld::CreateParticles() {" << endl;
    outfile << "#ifndef VECCORE_CUDA_DEVICE_COMPILATION\n";
    outfile << "  static bool fgCreateParticlesInitDone = false;\n";
    outfile << "#else\n";
    outfile << "  bool &fgCreateParticlesInitDone(fgCreateParticlesInitDoneDev);\n";
    outfile << "#endif\n";
    outfile << "  if (fgCreateParticlesInitDone) return;\n";
    outfile << "  fgCreateParticlesInitDone = true;\n";

    for (unsigned int i = 0; i < nfiles; ++i)
      outfile << "  CreateParticle" << setfill('0') << setw(4) << i << "();" << endl;
    outfile << "} // End of CreateParticle" << endl;
    outfile << "} // End of inline namespace" << endl;
    outfile << "} // End of geant namespace" << endl;
    outfile << "#if defined(__clang__) && !defined(__APPLE__)" << endl;
    outfile << "#pragma clang optimize on" << endl;
    outfile << "#endif" << endl;
    outfile.close();
  }
}
//________________________________________________________________________________________________
void ParticleOld::GetDecay(const std::string &line, int &dcount, Decay &decay)
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
void ParticleOld::NormDecay()
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
void ParticleOld::GetPart(const std::string &line, int &count, std::string &name, int &pdg, bool &matter, int &pcode,
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
std::ostream &operator<<(std::ostream &os, const ParticleOld::Decay &dec)
{
  os << "Type " << static_cast<int>(dec.fType) << " br " << dec.fBr << " products: ";
  for (unsigned int i = 0; i < dec.fDaughters.size(); ++i)
    os << " " << dec.fDaughters[i];

  /*   for(int i=0; i<dec.fDaughters.size(); ++i)
 os << " " << ParticleOld::Particles().at(dec.fDaughters[i]).Name(); */

  return os;
}
#endif
}
}
