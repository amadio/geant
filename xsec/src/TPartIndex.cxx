#include <TPartIndex.h>
#include <string>
#include <algorithm>
#include <TBuffer.h>
#include <iostream>

using std::transform;
using std::string;
using std::map;

#ifdef USE_ROOT
ClassImp(TPartIndex)
#endif

    const
    char *TPartIndex::fgPrName[FNPROC] = {"Transport",    "MultScatt",   "Ionisation", "Decay",      "inElastic",
                                          "Elastic",      "RestCapture", "Brehms",     "PairProd",   "Annihilation",
                                          "CoulombScatt", "Photoel",     "Compton",    "Conversion", "Capture",
                                          "Fission",      "Killer",      "Total"};
const short TPartIndex::fgPCode[FNPROC] = {1091, 2010, 2002, 6201, 4121, 4111, 4151, 2003, 2004,
                                           2005, 2001, 2012, 2013, 2014, 4131, 4141, 7403, 999};

const char *TPartIndex::fgEleSymbol[NELEM] = {
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na",  "Mg", "Al",  "Si", "P",   "S",  "Cl",
    "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni",  "Cu", "Zn",  "Ga", "Ge",  "As", "Se",
    "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",  "Pd", "Ag",  "Cd", "In",  "Sn", "Sb",
    "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",  "Eu", "Gd",  "Tb", "Dy",  "Ho", "Er",
    "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au",  "Hg", "Tl",  "Pb", "Bi",  "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm",  "Bk", "Cf",  "Es", "Fm",  "Md", "No",
    "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Uut", "Fl", "Uup", "Lv", "Uus", "Uuo"};
const char *TPartIndex::fgEleName[NELEM] = {
    "Hydrogen",    "Helium",    "Lithium",      "Beryllium",   "Boron",        "Carbon",        "Nitrogen",
    "Oxygen",      "Fluorine",  "Neon",         "Sodium",      "Magnesium",    "Aluminium",     "Silicon",
    "Phosphorus",  "Sulfur",    "Chlorine",     "Argon",       "Potassium",    "Calcium",       "Scandium",
    "Titanium",    "Vanadium",  "Chromium",     "Manganese",   "Iron",         "Cobalt",        "Nickel",
    "Copper",      "Zinc",      "Gallium",      "Germanium",   "Arsenic",      "Selenium",      "Bromine",
    "Krypton",     "Rubidium",  "Strontium",    "Yttrium",     "Zirconium",    "Niobium",       "Molybdenum",
    "Technetium",  "Ruthenium", "Rhodium",      "Palladium",   "Silver",       "Cadmium",       "Indium",
    "Tin",         "Antimony",  "Tellurium",    "Iodine",      "Xenon",        "Caesium",       "Barium",
    "Lanthanum",   "Cerium",    "Praseodymium", "Neodymium",   "Promethium",   "Samarium",      "Europium",
    "Gadolinium",  "Terbium",   "Dysprosium",   "Holmium",     "Erbium",       "Thulium",       "Ytterbium",
    "Lutetium",    "Hafnium",   "Tantalum",     "Tungsten",    "Rhenium",      "Osmium",        "Iridium",
    "Platinum",    "Gold",      "Mercury",      "Thallium",    "Lead",         "Bismuth",       "Polonium",
    "Astatine",    "Radon",     "Francium",     "Radium",      "Actinium",     "Thorium",       "Protactinium",
    "Uranium",     "Neptunium", "Plutonium",    "Americium",   "Curium",       "Berkelium",     "Californium",
    "Einsteinium", "Fermium",   "Mendelevium",  "Nobelium",    "Lawrencium",   "Rutherfordium", "Dubnium",
    "Seaborgium",  "Bohrium",   "Hassium",      "Meitnerium",  "Darmstadtium", "Roentgenium",   "Copernicium",
    "Ununtrium",   "Flerovium", "Ununpentium",  "Livermorium", "Ununseptium",  "Ununoctium"};

const float TPartIndex::fgWElem[NELEM] = {
    1.008,  4.0026, 6.94,   9.0122, 10.81,  12.011, 14.007, 15.999, 18.998, 20.180, 22.990, 24.305, 26.982, 28.085,
    30.974, 32.06,  35.45,  39.948, 39.098, 40.078, 44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693,
    63.546, 65.38,  69.723, 72.63,  74.922, 78.96,  79.904, 83.798, 85.468, 87.62,  88.906, 91.224, 92.906, 95.96,
    97.91,  101.07, 102.91, 106.42, 107.87, 112.41, 114.82, 118.71, 121.76, 127.60, 126.90, 131.29, 132.91, 137.33,
    138.91, 140.12, 140.91, 144.24, 144.91, 150.36, 151.96, 157.25, 158.93, 162.50, 164.93, 167.26, 168.93, 173.05,
    174.97, 178.49, 180.95, 183.84, 186.21, 190.23, 192.22, 195.08, 196.97, 200.59, 204.38, 207.2,  208.98, 208.98,
    209.99, 222.02, 223.02, 226.03, 227.03, 232.04, 231.04, 238.03, 237.05, 244.06, 243.06, 247.07, 247.07, 251.08,
    252.08, 257.10, 258.10, 259.10, 262.11, 265.12, 268.13, 271.13, 270,    277.15, 276.15, 281.16, 280.16, 285.17,
    284.18, 289.19, 288.19, 293,    294,    294};

TPartIndex *TPartIndex::fgPartIndex = 0;

//___________________________________________________________________
TPartIndex::TPartIndex()
    : fNPart(0), fPDG(0), fNpReac(0), fNpCharge(0), fNEbins(0), fEilDelta(0), fEGrid(0),
#ifdef USE_VECGEOM_NAVIGATOR
      fDBPdg(0),
#else
      fDBPdg(TDatabasePDG::Instance()),
#endif
      fPDGToGVMap() {
}

//___________________________________________________________________
TPartIndex::~TPartIndex() {
  delete[] fPDG;
  delete[] fEGrid;
  delete fDBPdg;
  fgPartIndex = 0;
}

//___________________________________________________________________
void TPartIndex::SetEnergyGrid(double emin, double emax, int nbins) {
  fNEbins = nbins;
  fEilDelta = (fNEbins - 1) / log(emax / emin);
  delete[] fEGrid;
  fEGrid = new double[fNEbins];
  double en = emin;
  double edelta = exp(1 / fEilDelta);
  for (int i = 0; i < fNEbins; ++i) {
    fEGrid[i] = en;
    en *= edelta;
  }
}

//___________________________________________________________________
int TPartIndex::ProcIndex(int proccode) const {
  short ip = fgNProc;
  while (ip--)
    if (fgPCode[ip] == proccode)
      break;
  return ip;
}

//___________________________________________________________________
const char *TPartIndex::ProcName(int proc) const {
  if (proc < 0 || proc >= fgNProc)
    return "Unknown";
  return fgPrName[proc];
}

//___________________________________________________________________
void TPartIndex::SetPartTable(const int *vpdg, int np) {
  fNPart = np;
  delete[] fPDG;
  fPDG = new int[fNPart];
  for (int i = 0; i < fNPart; ++i)
    fPDG[i] = vpdg[i];
}

//______________________________________________________________________________
int TPartIndex::PDG(const char *pname) const {
  int nr = fNPart;
#ifdef USE_VECGEOM_NAVIGATOR
  while (nr--)
    if (!strcmp(pname, fDBPdg->GetParticle(fPDG[nr]).Name()))
      return fPDG[nr];
#else
  while (nr--)
    if (!strcmp(pname, fDBPdg->GetParticle(fPDG[nr])->GetName()))
      return fPDG[nr];
#endif
  return -12345678;
}

//______________________________________________________________________________
int TPartIndex::PartIndex(int pdg) const {
  switch (pdg) {
  case 11:
    return fSpecGVIndices[0]; // e-
    break;
  case -11:
    return fSpecGVIndices[1]; // e+
    break;
  case 22:
    return fSpecGVIndices[2]; // gamma
    break;
  case 2212:
    return fSpecGVIndices[3]; // proton
    break;
  default: // look for in the map
    return fPDGToGVMap.find(pdg)->second;
  }
}

//______________________________________________________________________________
void TPartIndex::Print(const char *option) const {
  char line[120];
  string opt = option;
  transform(opt.begin(), opt.end(), opt.begin(), ::tolower);
  if (opt.find("particles") != string::npos) {
    printf("Available particles:\n");
    memset(line, 0, 120);
    for (int i = 0; i < fNPart; ++i) {
#ifdef USE_VECGEOM_NAVIGATOR
      const char *name = fDBPdg->GetParticle(fPDG[i]).Name();
#else
      const char *name = fDBPdg->GetParticle(fPDG[i])->GetName();
#endif
      if (strlen(line) + strlen(name) + 1 > 119) {
        printf("%s\n", line);
        memset(line, 0, 120);
      }
      strcat(line, name);
      strcat(line, " ");
    }
    if (strlen(line))
      printf("%s\n", line);
  }
  if (opt.find("reactions") != string::npos) {
    printf("Available reactions:\n");
    memset(line, 0, 120);
    strcat(line, "Total ");
    for (int i = 0; i < fgNProc; ++i) {
      if (strlen(line) + strlen(fgPrName[i]) + 1 > 119) {
        printf("%s\n", line);
        memset(line, 0, 120);
      }
      strcat(line, fgPrName[i]);
      strcat(line, " ");
    }
    if (strlen(line))
      printf("%s\n", line);
  }
  if (opt.find("version") != string::npos) {
    printf("Xsec database version %d.%d.%d\n", VersionMajor(), VersionMinor(), VersionSub());
  }
}

//______________________________________________________________________________
void TPartIndex::SetPDGToGVMap(std::map<int, int> &theMap) {
  fPDGToGVMap = theMap;
  fSpecGVIndices[0] = fPDGToGVMap.find(11)->second;   // e-
  fSpecGVIndices[1] = fPDGToGVMap.find(-11)->second;  // e+
  fSpecGVIndices[2] = fPDGToGVMap.find(22)->second;   // gamma
  fSpecGVIndices[3] = fPDGToGVMap.find(2212)->second; // proton
}

#ifdef USE_ROOT
//______________________________________________________________________________
void TPartIndex::Streamer(TBuffer &R__b) {
  // Stream an object of class TPartIndex.

  if (R__b.IsReading()) {
#ifdef USE_VECGEOM_NAVIGATOR
    delete fDBPdg;
    fDBPdg = 0;
#endif
    R__b.ReadClassBuffer(TPartIndex::Class(), this);
    fgPartIndex = this;

    Print("version");
    // create the particle reference vector
    fGVParticle.resize(fPDGToGVMap.size(), 0);
    for (map<int, int>::iterator p = fPDGToGVMap.begin(); p != fPDGToGVMap.end(); ++p) {
// std::cout << " gv index " << p->second << " corresponds to " << p->first << std::endl;
// create direct access vector with GeantV code
#ifdef USE_VECGEOM_NAVIGATOR
      const Particle_t *pp = &Particle_t::GetParticle(p->first);
      // set the code inside the particle too
      const_cast<Particle_t *>(pp)->SetCode(p->second);
      if (pp->Mass() >= 0)
        fGVParticle[p->second] = pp;
      else
        std::cout << __func__ << "::"
                  << " particle PDG " << p->first << " not found !" << std::endl;
#else
      const Particle_t *pp = TDatabasePDG::Instance()->GetParticle(p->first);
      if (!pp)
        std::cout << " Particle " << p->first << " does not exist " << std::endl;
      fGVParticle[p->second] = pp;
#ifdef SPECIAL_FCA_HACK
      //
      // this is a particularly ugly piece of code to write a text file (in fact stdout)
      // in the format of the root pdg_table.txt to add the particles defined by Geant4
      // to those defined by ROOT. Most of these particles are of little practical use
      // most probably, however this is mostly to avoid problems moving from GeantV to Geant4
      //
      if (p->second == 1) {
        int count = 521;
        // list of PDGs in Geant4 and not in ROOT
        // Some are only defined by the antiparticle, but it is enough to catch them all
        //
        int tupdg[231] = {
            -1000020040, -1000020030, -1000010030, -1000010020, -100012210, -100012110, -100002210, -100002110,
            -9000211,    -100325,     -100323,     -100321,     -100315,    -100313,    -100311,    -100213,
            -100211,     -53122,      -52214,      -52114,      -43122,     -42212,     -42124,     -42112,
            -41214,      -33324,      -33314,      -33122,      -32224,     -32214,     -32212,     -32124,
            -32114,      -32112,      -31214,      -31114,      -30323,     -30313,     -30213,     -23324,
            -23314,      -23224,      -23222,      -23214,      -23212,     -23126,     -23124,     -23122,
            -23114,      -23112,      -22224,      -22222,      -22214,     -22212,     -22124,     -22122,
            -22114,      -22112,      -21214,      -21212,      -21114,     -21112,     -13326,     -13324,
            -13316,      -13314,      -13226,      -13224,      -13222,     -13216,     -13214,     -13212,
            -13126,      -13124,      -13122,      -13116,      -13114,     -13112,     -12226,     -12224,
            -12222,      -12218,      -12216,      -12214,      -12212,     -12126,     -12122,     -12118,
            -12116,      -12114,      -12112,      -11216,      -11212,     -11116,     -11114,     -11112,
            -10325,      -10315,      -10215,      -3228,       -3226,      -3218,      -3216,      -3128,
            -3126,       -3124,       -3118,       -3116,       -2228,      -2226,      -2222,      -2218,
            -2216,       -2128,       -2126,       -2124,       -2122,      -2118,      -2116,      -1218,
            -1216,       -1214,       -1212,       -1118,       -1116,      -1112,      -327,       -317,
            -217,        117,         217,         227,         317,        327,        337,        1112,
            1116,        1118,        1212,        1214,        1216,       1218,       2116,       2118,
            2122,        2124,        2126,        2128,        2216,       2218,       2222,       2226,
            2228,        3116,        3118,        3124,        3126,       3128,       3216,       3218,
            3226,        3228,        10115,       10215,       10225,      10315,      10325,      10335,
            11112,       11114,       11116,       11212,       11216,      12112,      12114,      12116,
            12118,       12122,       12126,       12212,       12214,      12216,      12218,      12222,
            12224,       12226,       13112,       13114,       13116,      13122,      13124,      13126,
            13212,       13214,       13216,       13222,       13224,      13226,      13314,      13316,
            13324,       13326,       21112,       21114,       21212,      21214,      22112,      22114,
            22122,       22124,       22212,       22214,       22222,      30113,      30223,      100111,
            100113,      100221,      100223,      100331,      100333,     9000111,    9000221,    9010221,
            9020221,     9030221,     9030225,     9060225,     50000050,   50000052,   50000060};
        int updg[1000];
        int tot = 231;
        // Make sure that we add the particle and the antiparticle, extend the list if
        // necessary
        for (int l = 0; l < 231; ++l) {
          updg[l] = tupdg[l];
          bool found = false;
          for (int k = 0; k < 231; ++k)
            if (updg[l] == -tupdg[k]) {
              found = true;
              break;
            }
          if (!found)
            updg[tot++] = -tupdg[l];
          TParticlePDG *ap = TDatabasePDG::Instance()->GetParticle(-updg[l]);
          if (!ap) {
            // massage the name, we want antiparticles to have _bar at the end
            // as per the root convention and not anti_ at the beginnig as
            // Geant4
            TParticlePDG *p = TDatabasePDG::Instance()->GetParticle(updg[l]);
            string name(p->GetName());
            int index = name.find("Anti");
            if (index != std::string::npos)
              name.replace(index, 4, "");
            index = name.find("anti_");
            if (index != std::string::npos)
              name.replace(index, 5, "");
            if (updg[l] > 0)
              name += "_bar";
            // add the (anti)particle
            // Note that I discovered with horror that TParticlePDG has isospin, iso3 and
            // strangeness, but there is no way to set or to get them, we'll fix this later
            TDatabasePDG::Instance()->AddParticle(name.c_str(), p->GetTitle(), p->Mass(), p->Stable(), p->Width(),
                                                  -p->Charge(), p->ParticleClass(), -p->PdgCode(),
                                                  p->PdgCode() > 0 ? 1 : 0, p->TrackingCode());
          }
        }
        // Two ordering loops to make sure that the particle is always before
        // antiparticle as per pdg_table.txt conventions
        for (int l = 0; l < tot - 1; ++l)
          for (int k = l + 1; k < tot; ++k)
            if (updg[l] < updg[k]) {
              int u = updg[l];
              updg[l] = updg[k];
              updg[k] = u;
            }
        for (int l = 0; l < tot - 1; ++l)
          for (int k = l + 1; k < tot; ++k)
            if (abs(updg[l]) < abs(updg[k]) || (updg[l] == -updg[k] && updg[l] < updg[k])) {
              int u = updg[l];
              updg[l] = updg[k];
              updg[k] = u;
            }
        // Here actually write out the new lines
        // should really have been done on a file
        // but since I am not planning to use this code more than
        // once, I just copy-pasted from the screen
        for (int i = 0; i < tot; ++i) {
          count++;
          TParticlePDG *d = TDatabasePDG::Instance()->GetParticle(updg[i]);
          int ap = 1;
          string name(d->GetName());
          char fmt[100];
          int lpdg = log10(abs(updg[i])) + 1;
          if (updg[i] < 0) {
            ap = 0;
            lpdg += 1;
            for (int j = 0; j < tot; ++j)
              if (updg[i] == -updg[j]) {
                ap = j + 522;
                break;
              }
            if (ap == 0)
              exit(1);
            int index = name.find("Anti");
            if (index != std::string::npos)
              name.replace(index, 4, "");
            index = name.find("anti_");
            if (index != std::string::npos)
              name.replace(index, 5, "");
            index = name.find("_bar");
            if (index == std::string::npos)
              name += "_bar";
            int acode = 0;
            int len = std::max<int>(name.size(), 28 - lpdg - 1);
            sprintf(fmt, "%%5d %%-%d.%ds%%-%dd%%4d%%6d%%s", len, len, lpdg + 1);
            //		printf("%5d %-15.15s %-12d%4d%6d\n",
            printf(fmt, count, name.c_str(), updg[i], ap, acode, "\n");
          } else {
            int code = 100;
            string pclass = "Unknown";
            if (updg[i] > 1000000000)
              pclass = "ion";
            int iso = -100;
            if (d->Isospin())
              iso = d->Isospin() * 2 + 1;
            int s = -100;
            if (d->Strangeness())
              s = 2 * d->Strangeness() + 1;
            int flav = -1;
            int trkcode = -1;
            int ndec = 0;
            int len = 23 - lpdg - 1;
            len = std::max<int>(name.size(), 23 - lpdg - 1);
            sprintf(fmt, "%%5d %%-%d.%ds%%%dd%%3d%%4d %%-11s%%3d%%12.5e%%12.5e%%5d%%3d%%5d%%3d%%5d%%4d%%s", len, len,
                    lpdg + 1);
            //		printf("%s\n",fmt);
            // printf("%5d %-13.13s%10d%3d%4d %-10s%3d%12.5e%12.5e%5d%3d%5d%3d%5d%4d\n",
            printf(fmt, count, name.c_str(), updg[i], ap, code, pclass.c_str(), (int)d->Charge(), d->Mass(), d->Width(),
                   iso, (int)d->I3(), s, flav, trkcode, ndec, "\n");
          }
        }
      }
#endif
#endif
    }
  } else {
    R__b.WriteClassBuffer(TPartIndex::Class(), this);
  }
}
#endif

// semi empirical Bethe-Weizsacker mass formula based on the liquid drop model
// with coefficients determined by fitting to experimental mass data (AME2003)and
// taken from J.M. Pearson: Nuclear Mass Models (jina web /pearson_MassSchool)
// accuracy about +-few MeV
//______________________________________________________________________________
double TPartIndex::GetAprxNuclearMass(int Z, int A) {
  const double a_vol = -0.01569755; // volume    term coef. [GeV]
  const double a_surf = 0.01766269; // surface   term coef. [GeV]
  const double a_c = 0.0007070236;  // Coulomb   term coef. [GeV]
  const double a_sym = 0.026308165; // asymmetry term coef. [GeV]
  const double a_ss = -0.017003132; // surface-sym. term coef. [GeV]
  const double massp = 0.938272;    // mass of proton  [GeV]
  const double massn = 0.939565;    // mass of neutron [GeV]

  int N = A - Z;      // #neutrons
  double delta = 0.0; // if A is odd

  if ((N & 1) && (Z & 1)) // (N,Z) are odd
    delta = 0.0025;
  else if ((!(N & 1)) && (!(Z & 1))) // (N,Z) are even
    delta = -0.0025;

  double Ap13 = pow(A, 1 / 3.); // A^{1/3}
  double Ap23 = Ap13 * Ap13;    // A^{2/3}
  double dummy = (N - Z) / (1. * A);

  return Z * massp + N * massn +
         (a_vol * A + a_surf * Ap23 + a_c * Z * Z / Ap13 + (a_sym * A + a_ss * Ap23) * dummy * dummy + delta);
}
