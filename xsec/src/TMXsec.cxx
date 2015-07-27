#include "TMXsec.h"
#include "TPartIndex.h"
#include "TPDecay.h"
#include "TMath.h"
#include "TRandom.h"
#include "GeantTrack.h"
#include "GeantPropagator.h"
#include "GeantTaskData.h"
#ifdef USE_VECGEOM_NAVIGATOR
#include "volumes/Particle.h"
using vecgeom::Particle;
#endif
#include "base/Global.h"
using vecgeom::kAvogadro;

#include <algorithm>
using std::numeric_limits;
using std::max;

ClassImp(TMXsec)

    //____________________________________________________________________________
    TMXsec::TMXsec()
    : fNEbins(0), fNTotXL(0), fNCharge(0), fNRelXS(0), fEilDelta(0), fEGrid(0), fNElems(0), fElems(0), fTotXL(0),
      fRelXS(0), fDEdx(0), fMSangle(0), fMSansig(0), fMSlength(0), fMSlensig(0), fRatios(0), fRange(0),
      fInvRangeTable(0), fDecayTable(0) {
  fName[0] = '\0';
  fTitle[0] = '\0';
}

//____________________________________________________________________________
TMXsec::~TMXsec() {
  delete[] fElems;
  delete[] fTotXL;
  delete[] fRelXS;
  delete[] fMSangle;
  delete[] fMSansig;
  delete[] fMSlength;
  delete[] fMSlensig;
  delete[] fDEdx;
  delete[] fRatios;
  delete[] fRange;
  for (int i = 0; i < TPartIndex::I()->NPartCharge(); ++i)
    if (fInvRangeTable[i])
      delete[] fInvRangeTable[i];
}

//____________________________________________________________________________
TMXsec::TMXsec(const char *name, const char *title, const int z[], const int /*a*/[], const float w[], int nel,
               float dens, bool weight, const TPDecay *decaytable)
    : fNEbins(0), fNTotXL(0), fNCharge(0), fNRelXS(0), fEilDelta(TPartIndex::I()->EilDelta()),
      fEGrid(TPartIndex::I()->EGrid()), fNElems(0), fElems(0), fTotXL(0), fRelXS(0), fDEdx(0), fMSangle(0), fMSansig(0),
      fMSlength(0), fMSlensig(0), fRatios(0), fRange(0), fInvRangeTable(0), fDecayTable(0) {
  // Create a mixture material, we support only natural materials for the moment
  // so we ignore a (i.e. we consider it == 0)

  strncpy(fName, name, 31);
  strncpy(fTitle, title, 127);

  fDecayTable = decaytable;
  fNElems = nel;
  fElems = new TEXsec *[fNElems];
  memset(fElems, 0, fNElems * sizeof(TEXsec *));

  for (int i = 0; i < fNElems; ++i)
    if (z[i])
      fElems[i] = TEXsec::GetElement(z[i], 0);
    else if (fNElems > 1)
      Fatal("TMXsec", "Cannot have vacuum in mixtures");

  if (!z[0])
    return;

  fNEbins = TPartIndex::I()->NEbins();
  // double *ratios = new double[fNElems];
  fRatios = new double[fNElems];
  double *rdedx = new double[fNElems];
  double hnorm = 0;
  if (fNElems > 1) {
    for (int i = 0; i < fNElems; ++i) {
      // ratios[i] = w[i];
      // if(weight) ratios[i]/=TPartIndex::I()->WEle(z[i]);
      // hnorm+=ratios[i]*TPartIndex::I()->WEle(z[i]);
      fRatios[i] = w[i];
      if (weight)
        fRatios[i] /= TPartIndex::I()->WEle(z[i]);
      hnorm += fRatios[i] * TPartIndex::I()->WEle(z[i]);
    }
  } else {
    // ratios[0]=1;
    fRatios[0] = 1.0;
    hnorm = TPartIndex::I()->WEle(z[0]);
  }

  //   if(weight) printf("By weight: ");
  //   else       printf("By number: ");

  for (int i = 0; i < fNElems; ++i) {
    // rdedx[i] = ratios[i]*dens/fElems[i]->Dens();
    // ratios[i]*=kAvogadro*1e-24*dens/hnorm;
    if (fNElems > 1)
      rdedx[i] = fRatios[i] * TPartIndex::I()->WEle(z[i]) / hnorm; // mass fraction
    else
      rdedx[i] = 1.0;
    fRatios[i] *= kAvogadro * 1e-24 * dens / hnorm;
    //      printf("%d %f ",z[i],ratios[i]);
  }
  //   printf("\n");

  // Build table with total x-sections for all mate / parts

  int totindex = TPartIndex::I()->ProcIndex("Total");
  int npart = TPartIndex::I()->NPartReac();
  int ncharge = TPartIndex::I()->NPartCharge();
  // Layout part1 { en<1> { tot<1>, ... , tot<fNElems>}, .....en<nbins> {tot<1>, ..., tot<fNElems>}}

  if (fNElems > 1) {
    fNRelXS = npart * fNEbins * fNElems;
    fRelXS = new float[fNRelXS];
    memset(fRelXS, 0, fNRelXS * sizeof(float));
  }
  fNTotXL = npart * fNEbins;
  fTotXL = new float[fNTotXL];
  memset(fTotXL, 0, fNTotXL * sizeof(float));
  fNCharge = ncharge * fNEbins;
  fMSangle = new float[fNCharge];
  fMSansig = new float[fNCharge];
  fMSlength = new float[fNCharge];
  fMSlensig = new float[fNCharge];
  fDEdx = new float[fNCharge];
  memset(fMSangle, 0, fNCharge * sizeof(float));
  memset(fMSansig, 0, fNCharge * sizeof(float));
  memset(fMSlength, 0, fNCharge * sizeof(float));
  memset(fMSlensig, 0, fNCharge * sizeof(float));
  memset(fDEdx, 0, fNCharge * sizeof(float));

  for (int ip = 0; ip < npart; ++ip) {
    int ibase = ip * (fNEbins * fNElems);
    for (int ie = 0; ie < fNEbins; ++ie) {
      int ibin = ibase + ie * fNElems;
      if (fNElems > 1) {
        for (int iel = 0; iel < fNElems; ++iel) {
          fRelXS[ibin + iel] = fElems[iel]->XS(ip, totindex, fEGrid[ie]) * fRatios[iel];
          fTotXL[ip * fNEbins + ie] += fRelXS[ibin + iel];
        }
      } else {
        fTotXL[ip * fNEbins + ie] = fElems[0]->XS(ip, totindex, fEGrid[ie]) * fRatios[0];
      }
      if (fTotXL[ip * fNEbins + ie]) {
        fTotXL[ip * fNEbins + ie] = 1. / fTotXL[ip * fNEbins + ie];
        if (fNElems > 1)
          for (int iel = 0; iel < fNElems; ++iel)
            fRelXS[ibin + iel] *= fTotXL[ip * fNEbins + ie];
      } // else {
      //  fTotXL[ip*fNEbins+ie]=0.0;//numeric_limits<float>::max();
      // }
    }
  }

// add contribution from decay to the total mean free path of particles that
// have reactions other than decay; in case of particles that have only decay
// the decay mean free path will be computed on the fly
#ifdef USE_VECGEOM_NAVIGATOR
  if (Particle::GetParticle(11).Mass() <= 0)
    Particle::CreateParticles();
#endif
  for (int ip = 0; ip < npart; ++ip) {
    if (fDecayTable->HasDecay(ip)) {
      int pdgcode = TPartIndex::I()->PDG(ip);
#ifdef USE_VECGEOM_NAVIGATOR
      const Particle *const &partPDG = &Particle::GetParticle(pdgcode);
#else
      TParticlePDG *partPDG = TPartIndex::I()->DBPdg()->GetParticle(pdgcode);
#endif
      double mass = partPDG->Mass();                        // mass of the particle [GeV]
      double cTauPerMass = fDecayTable->GetCTauPerMass(ip); // c*tau/mass [cm/GeV]
      double lambda = 0.0;
      for (int ie = 0; ie < fNEbins; ++ie) {
        // decay mean free path [cm]: Ptotal*c*tau/mass
        lambda = sqrt(fEGrid[ie] * (fEGrid[ie] + 2.0 * mass)) * cTauPerMass;
        if (fTotXL[ip * fNEbins + ie] > 0.0)
          fTotXL[ip * fNEbins + ie] = fTotXL[ip * fNEbins + ie] * lambda / (fTotXL[ip * fNEbins + ie] + lambda);
        else
          fTotXL[ip * fNEbins + ie] = lambda;
      }
    } else {
      for (int ie = 0; ie < fNEbins; ++ie)
        if (fTotXL[ip * fNEbins + ie] == 0.0)
          fTotXL[ip * fNEbins + ie] = numeric_limits<float>::max();
    }
  }

  for (int ip = 0; ip < ncharge; ++ip) {
    float ang;
    float asig;
    float len;
    float lsig;
    for (int ie = 0; ie < fNEbins; ++ie) {
      for (int iel = 0; iel < fNElems; ++iel) {
        fElems[iel]->MS(ip, fEGrid[ie], ang, asig, len, lsig);
        fMSangle[ip * fNEbins + ie] += ang * rdedx[iel];
        fMSansig[ip * fNEbins + ie] += asig * rdedx[iel];
        fMSlength[ip * fNEbins + ie] += len * rdedx[iel];
        fMSlensig[ip * fNEbins + ie] += lsig * rdedx[iel];
        fDEdx[ip * fNEbins + ie] +=
            fElems[iel]->DEdx(ip, fEGrid[ie]) * rdedx[iel] * dens; // dE/dx [GeV/cm]; Bragg's rule
      }
    }
  }

  // build range table with simple integration of 1/dedx and the corresponding
  // inverse range table
  fRange = new float[fNCharge];
  memset(fRange, 0, fNCharge * sizeof(float));
  fInvRangeTable = new std::vector<std::pair<float, double>> *[fNCharge];
  memset(fInvRangeTable, 0, fNCharge * sizeof(std::vector<std::pair<float, double>> *));
  for (int ip = 0; ip < ncharge; ++ip) {
    std::vector<std::pair<float, double>> *aTable = new std::vector<std::pair<float, double>>();
    for (int ie = 0; ie < fNEbins; ++ie) {
      if (ie == 0) {
        // assume linear dependece
        fRange[ip * fNEbins + ie] = 2.0 * fEGrid[0] / fDEdx[ip * fNEbins + 0];
      } else {
        double xx =
            0.5 * (1.0 / fDEdx[ip * fNEbins + ie - 1] + 1.0 / fDEdx[ip * fNEbins + ie]) * (fEGrid[ie] - fEGrid[ie - 1]);
        fRange[ip * fNEbins + ie] = fRange[ip * fNEbins + ie - 1] + xx; // in [cm]
      }
      // insert into particle inverse range table
      aTable->push_back(std::make_pair(fRange[ip * fNEbins + ie], fEGrid[ie]));
    }
    // insert into inverse range table and sort by the first element of pairs
    std::sort(aTable->begin(), aTable->end());
    fInvRangeTable[ip] = aTable;
  }

  // normalization of fRatios[] to 1
  if (fNElems == 1)
    fRatios[0] = 1.0;
  else {
    for (int i = 1; i < fNElems; ++i)
      fRatios[i] += fRatios[i - 1];

    for (int i = 0; i < fNElems; ++i)
      fRatios[i] /= fRatios[fNElems - 1];
  }

  // cleaning up
  // delete [] ratios;
  delete[] rdedx;
}

/*
//____________________________________________________________________________
bool TMXsec::Prune() {
   // Prune elements
   for(int iel=0; iel<TEXsec::NLdElems(); ++iel) TEXsec::Element(iel)->Prune();
   return kTRUE;
}
*/

//____________________________________________________________________________
float TMXsec::Xlength(int part, float en, double ptot) {
  if (part < TPartIndex::I()->NPartReac()) {
    en = en <= fEGrid[fNEbins - 1] ? en : fEGrid[fNEbins - 1] * 0.999;
    en = max<double>(en,fEGrid[0]);
    int ibin = log(en / fEGrid[0]) * fEilDelta;
    ibin = ibin < fNEbins - 1 ? ibin : fNEbins - 2;
    //     double en1 = fEmin*exp(fElDelta*ibin);
    //     double en2 = en1*fEDelta;
    double en1 = fEGrid[ibin];
    double en2 = fEGrid[ibin + 1];
    //    if(en1>en || en2<en) {
    //      Error("Xlength","Wrong bin %d in interpolation: should be %f < %f < %f\n",
    //            ibin, en1, en, en2);
    //      return numeric_limits<float>::max();
    //    }
    double xrat = (en2 - en) / (en2 - en1);
    double x = xrat * fTotXL[part * fNEbins + ibin] + (1 - xrat) * fTotXL[part * fNEbins + ibin + 1];
    return x;
  } else if (fDecayTable->HasDecay(part)) {
    return ptot * fDecayTable->GetCTauPerMass(part); // Ptot*c*tau/mass [cm]
  } else {
    return numeric_limits<float>::max();
  }
}

/*
//____________________________________________________________________________
bool TMXsec::Xlength_v(int npart, const int part[], const float en[], double lam[])
{
  double ene;
  for(int ip=0; ip<npart; ++ip) {
    if(part[ip]>=TPartIndex::I()->NPartReac() || !fTotXL)
      lam[ip]=numeric_limits<float>::max();
    else {
      ene=en[ip]<=fEGrid[fNEbins-1]?en[ip]:fEGrid[fNEbins-1]*0.999;
      ene=max<double>(en[ip],fEGrid[0]);
      int ibin = log(ene/fEGrid[0])*fEilDelta;
      ibin = ibin<fNEbins-1?ibin:fNEbins-2;
      //     double en1 = fEmin*exp(fElDelta*ibin);
      //     double en2 = en1*fEDelta;
      double en1 = fEGrid[ibin];
      double en2 = fEGrid[ibin+1];
      if(en1>ene || en2<ene) {
        Error("Xlength","Wrong bin %d in interpolation: should be %f < %f < %f\n",
              ibin, en1, ene, en2);
        lam[ip] = numeric_limits<float>::max();
      }
      double xrat = (en2-ene)/(en2-en1);
      lam[ip] = xrat*fTotXL[part[ip]*fNEbins+ibin]+(1-xrat)*fTotXL[part[ip]*fNEbins+ibin+1];
    }
  }
  return kTRUE;
}
*/

//______________________________________________________________________________
void TMXsec::ProposeStep(int ntracks, GeantTrack_v &tracks, GeantTaskData *td) {
  // Propose step for the first ntracks in the input vector of tracks and write to
  // tracks.fPstepV[]

  // tid-based rng: need $\{ R_i\}_i^{ntracks} \in \mathcal{U} \in [0,1]$
  double *rndArray = td->fDblArray;
  td->fRndm->RndmArray(ntracks, rndArray);

  for (int i = 0; i < ntracks; ++i) {
    int ipart = tracks.fGVcodeV[i];                   // GV particle index/code
    double energy = tracks.fEV[i] - tracks.fMassV[i]; // $E_{kin}$
    energy = energy <= fEGrid[fNEbins - 1] ? energy : fEGrid[fNEbins - 1] * 0.999;
    energy = max<double>(energy,fEGrid[0]);

    // continous step limit if any
    tracks.fEindexV[i] = -1; // Flag continous step limit
    double cx = numeric_limits<double>::max();
    if (ipart < TPartIndex::I()->NPartCharge()) {
      double range = Range(ipart, energy);
      cx = range;
      if (cx < 0.) {
        cx = numeric_limits<double>::max();
      } else {
        const double dRR = .2;
        const double finRange = .1;
        if (cx > finRange)
          cx = cx * dRR + finRange * ((1. - dRR) * (2.0 - finRange / cx));
      }
    }
    tracks.fPstepV[i] = cx; // set it to the cont. step limit and update later

    // discrete step limit
    if (ipart < TPartIndex::I()->NPartReac()) {       // can have different reactions + decay
      int ibin = log(energy / fEGrid[0]) * fEilDelta; // energy bin index
      ibin = ibin < fNEbins - 1 ? ibin : fNEbins - 2;

      double en1 = fEGrid[ibin];
      double en2 = fEGrid[ibin + 1];

      double xrat = (en2 - energy) / (en2 - en1);
      // get interpolated -(total mean free path)
      double x = xrat * fTotXL[ipart * fNEbins + ibin] + (1 - xrat) * fTotXL[ipart * fNEbins + ibin + 1];
      x *= -1. * log(rndArray[i]);
      if (x < cx) {
        tracks.fPstepV[i] = x;
        tracks.fEindexV[i] = 1000; // Flag NOT continous step limit
      }
    } else if (fDecayTable->HasDecay(ipart)) {                       // it has only decay
      double x = tracks.fPV[i] * fDecayTable->GetCTauPerMass(ipart); // Ptot*c*tau/mass [cm]
      x = -1. * x * log(rndArray[i]);
      if (x < cx) {
        tracks.fPstepV[i] = x;
        tracks.fEindexV[i] = 1000; // Flag NOT continous step limit
      }
    }
  }
}

//______________________________________________________________________________
void TMXsec::ProposeStepSingle(int i, GeantTrack_v &tracks, GeantTaskData *td) {
  // Propose step for a single track in the input vector of tracks and write to
  // tracks.fPstepV[]

  // tid-based rng: need $\{ R_i\}_i^{ntracks} \in \mathcal{U} \in [0,1]$
  double *rndArray = td->fDblArray;
  td->fRndm->RndmArray(1, rndArray);

  int ipart = tracks.fGVcodeV[i];                   // GV particle index/code
  double energy = tracks.fEV[i] - tracks.fMassV[i]; // $E_{kin}$
  energy = energy <= fEGrid[fNEbins - 1] ? energy : fEGrid[fNEbins - 1] * 0.999;
  energy = max<double>(energy,fEGrid[0]);

  // continous step limit if any
  tracks.fEindexV[i] = -1; // Flag continous step limit
  double cx = numeric_limits<double>::max();
  if (ipart < TPartIndex::I()->NPartCharge()) {
    double range = Range(ipart, energy);
    cx = range;
    if (cx < 0.) {
      cx = numeric_limits<double>::max();
    } else {
      const double dRR = .2;
      const double finRange = .1;
      if (cx > finRange)
        cx = cx * dRR + finRange * ((1. - dRR) * (2.0 - finRange / cx));
    }
  }
  tracks.fPstepV[i] = cx; // set it to the cont. step limit and update later

  // discrete step limit
  if (ipart < TPartIndex::I()->NPartReac()) {       // can have different reactions + decay
    int ibin = log(energy / fEGrid[0]) * fEilDelta; // energy bin index
    ibin = ibin < fNEbins - 1 ? ibin : fNEbins - 2;

    double en1 = fEGrid[ibin];
    double en2 = fEGrid[ibin + 1];

    double xrat = (en2 - energy) / (en2 - en1);
    // get interpolated -(total mean free path)
    double x = xrat * fTotXL[ipart * fNEbins + ibin] + (1 - xrat) * fTotXL[ipart * fNEbins + ibin + 1];
    x *= -1. * log(rndArray[0]);
    if (x < cx) {
      tracks.fPstepV[i] = x;
      tracks.fEindexV[i] = 1000; // Flag NOT continous step limit
    }
  } else if (fDecayTable->HasDecay(ipart)) {                       // it has only decay
    double x = tracks.fPV[i] * fDecayTable->GetCTauPerMass(ipart); // Ptot*c*tau/mass [cm]
    x = -1. * x * log(rndArray[0]);
    if (x < cx) {
      tracks.fPstepV[i] = x;
      tracks.fEindexV[i] = 1000; // Flag NOT continous step limit
    }
  }
}

// get MS angles : NOT USED CURRENTLY (MS is not active)
//____________________________________________________________________________
float TMXsec::MS(int ipart, float energy) {
  int ibin;

  if (ipart >= TPartIndex::I()->NPartCharge())
    return 0.;

  double adj_energy = energy <= fEGrid[fNEbins - 1] ? energy : fEGrid[fNEbins - 1] * 0.999;
  adj_energy = max<double>(adj_energy,fEGrid[0]);

  ibin = log(adj_energy / fEGrid[0]) * fEilDelta;
  ibin = ibin < fNEbins - 1 ? ibin : fNEbins - 2;

  double en1 = fEGrid[ibin];
  double en2 = fEGrid[ibin + 1];

  if (en1 > adj_energy || en2 < adj_energy) {
    Error("MS", "Wrong bin %d in interpolation: should be %g < %g < %g\n", ibin, en1, adj_energy, en2);
    return 0.;
  }

  double xrat = (en2 - adj_energy) / (en2 - en1);
  return xrat * fMSangle[ipart * fNEbins + ibin] + (1 - xrat) * fMSangle[ipart * fNEbins + ibin + 1];
}

//____________________________________________________________________________
float TMXsec::DEdx(int part, float en) {
  if (part >= TPartIndex::I()->NPartCharge() || !fDEdx)
    return 0;
  else {
    en = en <= fEGrid[fNEbins - 1] ? en : fEGrid[fNEbins - 1] * 0.999;
    en = max<double>(en,fEGrid[0]);
    int ibin = log(en / fEGrid[0]) * fEilDelta;
    ibin = ibin < fNEbins - 1 ? ibin : fNEbins - 2;
    //     double en1 = fEmin*exp(fElDelta*ibin);
    //     double en2 = en1*fEDelta;
    double en1 = fEGrid[ibin];
    double en2 = fEGrid[ibin + 1];
    //    if(en1>en || en2<en) {
    //      Error("DEdx","Wrong bin %d in interpolation: should be %f < %f < %f\n",
    //            ibin, en1, en, en2);
    //      return numeric_limits<float>::max();
    //    }

    double xrat = (en2 - en) / (en2 - en1);
    return xrat * fDEdx[part * fNEbins + ibin] + (1 - xrat) * fDEdx[part * fNEbins + ibin + 1];
  }
}

/*
//____________________________________________________________________________
bool TMXsec::DEdx_v(int npart, const int part[], const float en[], float de[]) {
  double ene;
  if(!fDEdx) {
    memset(de,0,npart*sizeof(float));
    return kFALSE;
  }
  for(int ip=0; ip<npart; ++ip) {
    if(part[ip]>=TPartIndex::I()->NPartCharge())
      de[ip]=0;
    else {
      ene=en[ip]<=fEGrid[fNEbins-1]?en[ip]:fEGrid[fNEbins-1]*0.999;
      ene=max<double>(en[ip],fEGrid[0]);
      int ibin = log(ene/fEGrid[0])*fEilDelta;
      ibin = ibin<fNEbins-1?ibin:fNEbins-2;
      //     double en1 = fEmin*exp(fElDelta*ibin);
      //     double en2 = en1*fEDelta;
      double en1 = fEGrid[ibin];
      double en2 = fEGrid[ibin+1];
      if(en1>ene || en2<ene) {
        Error("DEdx","Wrong bin %d in interpolation: should be %f < %f < %f\n",
              ibin, en1, ene, en2);
        de[ip]=numeric_limits<float>::max();
      }
      double xrat = (en2-ene)/(en2-en1);
      de[ip] = xrat*fDEdx[part[ip]*fNEbins+ibin]+(1-xrat)*fDEdx[part[ip]*fNEbins+ibin+1];
    }
  }
  return kTRUE;
}
*/

//____________________________________________________________________________
float TMXsec::Range(int part, float en) {
  if (part >= TPartIndex::I()->NPartCharge() || !fRange)
    return -1.0;
  else {
    en = en <= fEGrid[fNEbins - 1] ? en : fEGrid[fNEbins - 1] * 0.999;
    en = max<double>(en,fEGrid[0]);
    int ibin = log(en / fEGrid[0]) * fEilDelta;
    ibin = ibin < fNEbins - 1 ? ibin : fNEbins - 2;

    double en1 = fEGrid[ibin];
    double en2 = fEGrid[ibin + 1];
    //    if(en1>en || en2<en) {
    //      Error("Range","Wrong bin %d in interpolation: should be %f < %f < %f\n",
    //            ibin, en1, en, en2);
    //      return -1.0;
    //    }

    double xrat = (en2 - en) / (en2 - en1);
    return xrat * fRange[part * fNEbins + ibin] + (1 - xrat) * fRange[part * fNEbins + ibin + 1];
  }
}

//____________________________________________________________________________
double TMXsec::InvRange(int part, float step) {
  if (part >= TPartIndex::I()->NPartCharge() || !fRange)
    return 0.;

  double minr = fRange[part * fNEbins]; // min available range i.e. for E_0
  if (step >= minr) {
    /*    int i=0;
        while(fRange[part*fNEbins+i]<step) ++i;
        double r1 = fRange[part*fNEbins+i-1];
        double r2 = fRange[part*fNEbins+i];
        double xrat = (r2-step)/(r2-r1);
        return xrat*fEGrid[i-1]+(1-xrat)*fEGrid[i];
    */
    if (step >= fInvRangeTable[part]->back().first)
      return fInvRangeTable[part]->back().second;
    // if(x<=(*table)[0].first) return (*table)[0].second;

    std::vector<std::pair<float, double>>::iterator itl, ith;
    ith = std::lower_bound(fInvRangeTable[part]->begin(), fInvRangeTable[part]->end(), std::make_pair(step, 0.0));
    itl = ith;
    --itl;
    double rat = (ith->first - step) / (ith->first - itl->first);
    return rat * itl->second + (1 - rat) * ith->second;

  } else {
    double x = step / minr;
    return fEGrid[0] * x * x;
  }
}

// Compute along step energy loss for charged particles using linear loss aprx.
//____________________________________________________________________________
void TMXsec::Eloss(int ntracks, GeantTrack_v &tracks) {
  // -should be called only for charged particles (first fNPartCharge particle
  // in TPartIndex::fPDG[]); the case ipartindex>=fNPartCharge is handled now in the if
  // Compute energy loss for the first ntracks in the input vector and update
  // tracks.fEV, tracks.fPV and tracks.EdepV. If the particle is stopped set the
  // necessary at-rest process type if it has any.

  double energyLimit = GeantPropagator::Instance()->fEmin;
  for (int i = 0; i < ntracks; ++i) {
    int ipart = tracks.fGVcodeV[i]; // GV particle index/code
    tracks.fProcessV[i] = -1;       // init process index to -1 i.e. no process
    double dedx = 0.0;

    // just a check; can be removed if it is ensured before calling
    if (ipart >= TPartIndex::I()->NPartCharge()) // there is no energy loss nothing change
      continue;

    double energy = tracks.fEV[i] - tracks.fMassV[i]; // $E_{kin}$ [GeV]
    double range = Range(ipart, energy);
    if (tracks.fStepV[i] > range) {
      // Particle will stop
      tracks.fEdepV[i] += energy;       // put Ekin to edepo
      tracks.fEV[i] = tracks.fMassV[i]; // set Etotal = Mass i.e. Ekin = 0
      tracks.fPV[i] = 0.;               // set Ptotal = 0
      tracks.fStatusV[i] = kKilled;     // set status to killed
      // stopped: set proc. indx = -2 to inidicate that
      tracks.fProcessV[i] = -2; // possible at-rest process need to be called
      continue;
    }

    // get dedx value
    if (energy <= fEGrid[0])
      dedx = fDEdx[ipart * fNEbins]; // protections
    else if (energy >= fEGrid[fNEbins - 1])
      dedx = fDEdx[ipart * fNEbins + fNEbins - 1];
    else {                                            // regular case
      int ibin = log(energy / fEGrid[0]) * fEilDelta; // energy bin index
      ibin = ibin < fNEbins - 1 ? ibin : fNEbins - 2;
      double en1 = fEGrid[ibin];
      double en2 = fEGrid[ibin + 1];
      double xrat = (en2 - energy) / (en2 - en1);
      dedx = xrat * fDEdx[ipart * fNEbins + ibin] + (1 - xrat) * fDEdx[ipart * fNEbins + ibin + 1];
    }
    // Update energy and momentum
    double gammaold = tracks.Gamma(i);
    double bgold = sqrt((gammaold - 1) * (gammaold + 1));

    double edepo = tracks.fStepV[i] * dedx; // compute energy loss using linera loss aprx.
    if (edepo > 0.01 * energy)              // long step: eloss > 1% of initial energy
      edepo = energy - InvRange(ipart, range - tracks.fStepV[i]);

    double newEkin = energy - edepo; // new kinetic energy
    if (newEkin < energyLimit) {     // new Kinetic energy below tracking cut
      // Particle energy below threshold
      tracks.fEdepV[i] += energy;       // put Ekin to edepo
      tracks.fEV[i] = tracks.fMassV[i]; // set Etotal = Mass i.e. Ekin = 0
      tracks.fPV[i] = 0.;               // set Ptotal = 0
      tracks.fStatusV[i] = kKilled;     // set status to killed
                                        // check if the particle stopped or just went below the tracking cut
      if (newEkin <= 0.0)               // stopped: set proc. indx = -2 to inidicate that
        tracks.fProcessV[i] = -2;       // possible at-rest process need to be called
      // else : i.e. 0 < newEkin < energyLimit just kill the track cause tracking cut
    } else { // track further: update energy, momentum, process etc.
      // tracks.fProcessV[i] = TPartIndex::I()->ProcIndex("Ionisation");
      tracks.fProcessV[i] = 2;   // Ionisation
      tracks.fEdepV[i] += edepo; // add energy deposit
      tracks.fEV[i] -= edepo;    // update total energy
      double gammanew = tracks.Gamma(i);
      double bgnew = sqrt((gammanew - 1) * (gammanew + 1));
      double pnorm = bgnew / bgold;
      tracks.fPV[i] *= pnorm; // update total momentum
    }
  }
}

// Compute along step energy loss for charged particles using linear loss aprx.
//____________________________________________________________________________
void TMXsec::ElossSingle(int i, GeantTrack_v &tracks) {
  // -should be called only for charged particles (first fNPartCharge particle
  // in TPartIndex::fPDG[]); the case ipartindex>=fNPartCharge is handled now in the if
  // Compute energy loss for the first ntracks in the input vector and update
  // tracks.fEV, tracks.fPV and tracks.EdepV. If the particle is stopped set the
  // necessary at-rest process type if it has any.

  double energyLimit = GeantPropagator::Instance()->fEmin;
  int ipart = tracks.fGVcodeV[i]; // GV particle index/code
  tracks.fProcessV[i] = -1;       // init process index to -1 i.e. no process
  double dedx = 0.0;

  // just a check; can be removed if it is ensured before calling
  if (ipart >= TPartIndex::I()->NPartCharge()) // there is no energy loss nothing change
    return;

  double energy = tracks.fEV[i] - tracks.fMassV[i]; // $E_{kin}$ [GeV]
  double range = Range(ipart, energy);
  if (tracks.fStepV[i] > range) {
    // Particle will stop
    tracks.fEdepV[i] += energy;       // put Ekin to edepo
    tracks.fEV[i] = tracks.fMassV[i]; // set Etotal = Mass i.e. Ekin = 0
    tracks.fPV[i] = 0.;               // set Ptotal = 0
    tracks.fStatusV[i] = kKilled;     // set status to killed
    // stopped: set proc. indx = -2 to inidicate that
    tracks.fProcessV[i] = -2; // possible at-rest process need to be called
    return;
  }

  // get dedx value
  if (energy <= fEGrid[0])
    dedx = fDEdx[ipart * fNEbins]; // protections
  else if (energy >= fEGrid[fNEbins - 1])
    dedx = fDEdx[ipart * fNEbins + fNEbins - 1];
  else {                                            // regular case
    int ibin = log(energy / fEGrid[0]) * fEilDelta; // energy bin index
    ibin = ibin < fNEbins - 1 ? ibin : fNEbins - 2;
    double en1 = fEGrid[ibin];
    double en2 = fEGrid[ibin + 1];
    double xrat = (en2 - energy) / (en2 - en1);
    dedx = xrat * fDEdx[ipart * fNEbins + ibin] + (1 - xrat) * fDEdx[ipart * fNEbins + ibin + 1];
  }
  // Update energy and momentum
  double gammaold = tracks.Gamma(i);
  double bgold = sqrt((gammaold - 1) * (gammaold + 1));

  double edepo = tracks.fStepV[i] * dedx; // compute energy loss using linera loss aprx.
  if (edepo > 0.01 * energy)              // long step: eloss > 1% of initial energy
    edepo = energy - InvRange(ipart, range - tracks.fStepV[i]);

  double newEkin = energy - edepo; // new kinetic energy
  if (newEkin < energyLimit) {     // new Kinetic energy below tracking cut
    // Particle energy below threshold
    tracks.fEdepV[i] += energy;       // put Ekin to edepo
    tracks.fEV[i] = tracks.fMassV[i]; // set Etotal = Mass i.e. Ekin = 0
    tracks.fPV[i] = 0.;               // set Ptotal = 0
    tracks.fStatusV[i] = kKilled;     // set status to killed
                                      // check if the particle stopped or just went below the tracking cut
    if (newEkin <= 0.0)               // stopped: set proc. indx = -2 to inidicate that
      tracks.fProcessV[i] = -2;       // possible at-rest process need to be called
    // else : i.e. 0 < newEkin < energyLimit just kill the track cause tracking cut
  } else { // track further: update energy, momentum, process etc.
    // tracks.fProcessV[i] = TPartIndex::I()->ProcIndex("Ionisation");
    tracks.fProcessV[i] = 2;   // Ionisation
    tracks.fEdepV[i] += edepo; // add energy deposit
    tracks.fEV[i] -= edepo;    // update total energy
    double gammanew = tracks.Gamma(i);
    double bgnew = sqrt((gammanew - 1) * (gammanew + 1));
    double pnorm = bgnew / bgold;
    tracks.fPV[i] *= pnorm; // update total momentum
  }
}

//____________________________________________________________________________
TEXsec *TMXsec::SampleInt(int part, double en, int &reac, double ptot) {
  if (part < TPartIndex::I()->NPartReac()) {
    // first sample if deacy or something else if the particle can have decay
    bool isDecay = kFALSE;
    if (fDecayTable->HasDecay(part)) {
      double lambdaDecay = ptot * fDecayTable->GetCTauPerMass(part);
      double lambdaTotal = Xlength(part, en, ptot); // ptotal is not used now there
      // $P(decay) =\lambda_{Total}/\lambda_{decay}$
      if (gRandom->Rndm() < lambdaTotal / lambdaDecay)
        isDecay = kTRUE;
    }
    if (isDecay) {
      reac = 3; // decay
      return 0; // on nothing
    } else {
      en = en <= fEGrid[fNEbins - 1] ? en : fEGrid[fNEbins - 1] * 0.999;
      en = max<double>(en,fEGrid[0]);
      int ibin = log(en / fEGrid[0]) * fEilDelta;
      ibin = ibin < fNEbins - 1 ? ibin : fNEbins - 2;
      //      double en1 = fEmin*exp(fElDelta*ibin);
      //    double en2 = en1*fEDelta;
      double en1 = fEGrid[ibin];
      double en2 = fEGrid[ibin + 1];
      //        if(en1>en || en2<en) {
      //	   Error("SampleInt","Wrong bin %d in interpolation: should be %f < %f < %f\n",
      //	         ibin, en1, en, en2);
      //	   reac=-1;
      //	   return 0;
      //        }
      int iel = -1;
      if (fNElems == 1) {
        iel = 0;
      } else {
        double xrat = (en2 - en) / (en2 - en1);
        double xnorm = 1.;
        while (iel < 0) {
          double ran = xnorm * gRandom->Rndm();
          double xsum = 0;
          for (int i = 0; i < fNElems; ++i) {
            xsum += xrat * fRelXS[ibin * fNElems + i] + (1 - xrat) * fRelXS[(ibin + 1) * fNElems + i];
            if (ran <= xsum) {
              iel = i;
              break;
            }
          }
          xnorm = xsum;
        }
      }
      reac = fElems[iel]->SampleReac(part, en);
      return fElems[iel];
    }
  } else if (fDecayTable->HasDecay(part)) {
    reac = 3; // decay
    return 0; // on nothing
  } else {
    reac = -1; // nothing
    return 0;  // on nothing
  }
}

//____________________________________________________________________________
void TMXsec::SampleInt(int ntracks, GeantTrack_v &tracksin, GeantTaskData *td) {
  int nParticleWithReaction = TPartIndex::I()->NPartReac();

  // tid-based rng
  double *rndArray = td->fDblArray;
  td->fRndm->RndmArray(2 * ntracks, rndArray);

  for (int t = 0; t < ntracks; ++t) {
    /*
         if(tracksin.fEindexV[t] < 0){ // continous step limited this step
             tracksin.fProcessV[t] = -1; // nothing
             tracksin.fEindexV[t]  = -1; // on nothing
           continue;
         }
    */
    int ipart = tracksin.fGVcodeV[t];
    // if the particle can have both decay and something else
    if (ipart < nParticleWithReaction) {
      bool isDecay = kFALSE;
      double energy = tracksin.fEV[t] - tracksin.fMassV[t]; // Ekin [GeV]
      if (fDecayTable->HasDecay(ipart)) {
        double ptotal = tracksin.fPV[t]; // Ptotal [GeV]
        double lambdaDecay = ptotal * fDecayTable->GetCTauPerMass(ipart);
        // ptot is not used now in Xlength(ipart,energy,ptot) since ipart<nParticleWithReaction
        double lambdaTotal = Xlength(ipart, energy, ptotal);
        // $P(decay) =\lambda_{Total}/\lambda_{decay}$
        if (rndArray[t] < lambdaTotal / lambdaDecay)
          isDecay = kTRUE;
      }
      if (isDecay) {
        tracksin.fProcessV[t] = 3; // decay
        tracksin.fEindexV[t] = -1; // on nothing
      } else {
        // not decay but something else; sample what else on what target
        energy = energy <= fEGrid[fNEbins - 1] ? energy : fEGrid[fNEbins - 1] * 0.999;
        energy = max<double>(energy,fEGrid[0]);

        int ibin = log(energy / fEGrid[0]) * fEilDelta;
        ibin = ibin < fNEbins - 1 ? ibin : fNEbins - 2;

        double en1 = fEGrid[ibin];
        double en2 = fEGrid[ibin + 1];
        // 1.
        // this material/mixture (TMXsec object) is composed of fNElems elements
        int iel = -1; // index of elements in TEXsec** fElems  ([fNElems]) of this
        if (fNElems == 1) {
          iel = 0;
        } else { // sampling one element based on the elemntal realtive X-secs that
          // have been normalized in CTR at the Ebins!; energy interp. is
          // included now-> while loop is to avoid problems from interp.
          double xrat = (en2 - energy) / (en2 - en1);
          double xnorm = 1.;
          while (iel < 0) {
            double ran = xnorm * td->fRndm->Rndm();
            double xsum = 0;
            for (int i = 0; i < fNElems; ++i) { // simple sampling from discrete p.
              xsum += xrat * fRelXS[ibin * fNElems + i] + (1 - xrat) * fRelXS[(ibin + 1) * fNElems + i];
              if (ran <= xsum) {
                iel = i;
                break;
              }
            }
            xnorm = xsum;
          }
        }
        // at this point the index of the element is sampled:= iel
        // the corresponding TEXsec* is fElems[iel]

        // sample the reaction by using the TEXsec* that corresponds to the iel-th
        // element i.e. fElems[iel] for the current particle (with particle index of
        // ipart) at the current energy; will retrun with the reaction index;
        tracksin.fProcessV[t] = fElems[iel]->SampleReac(ipart, energy, rndArray[ntracks + t]);
        tracksin.fEindexV[t] = fElems[iel]->Index(); // index of the selected element in TTabPhysMrg::fElemXsec[]
      }
    } else if (fDecayTable->HasDecay(ipart)) {
      // only decay can happen because ipart>nParticleWithReaction
      tracksin.fProcessV[t] = 3;  // decay
      tracksin.fEindexV[t] = -1;  // on nothing
    } else {                      // nothing happens
      tracksin.fProcessV[t] = -1; // nothing
      tracksin.fEindexV[t] = -1;  // on nothing
    }
  }
}

//____________________________________________________________________________
void TMXsec::SampleSingleInt(int t, GeantTrack_v &tracksin, GeantTaskData *td) {
  // Sample the interaction for a single particle at a time.
  int nParticleWithReaction = TPartIndex::I()->NPartReac();

  // tid-based rng
  double *rndArray = td->fDblArray;
  td->fRndm->RndmArray(2, rndArray);
  /*
     if(tracksin.fEindexV[t] < 0){ // continous step limited this step
         tracksin.fProcessV[t] = -1; // nothing
         tracksin.fEindexV[t]  = -1; // on nothing
       return;
     }
  */
  int ipart = tracksin.fGVcodeV[t];
  // if the particle can have both decay and something else
  if (ipart < nParticleWithReaction) {
    bool isDecay = kFALSE;
    double energy = tracksin.fEV[t] - tracksin.fMassV[t]; // Ekin [GeV]
    if (fDecayTable->HasDecay(ipart)) {
      double ptotal = tracksin.fPV[t]; // Ptotal [GeV]
      double lambdaDecay = ptotal * fDecayTable->GetCTauPerMass(ipart);
      // ptot is not used now in Xlength(ipart,energy,ptot) since ipart<nParticleWithReaction
      double lambdaTotal = Xlength(ipart, energy, ptotal);
      // $P(decay) =\lambda_{Total}/\lambda_{decay}$
      if (rndArray[0] < lambdaTotal / lambdaDecay)
        isDecay = kTRUE;
    }
    if (isDecay) {
      tracksin.fProcessV[t] = 3; // decay
      tracksin.fEindexV[t] = -1; // on nothing
    } else {
      // not decay but something else; sample what else on what target
      energy = energy <= fEGrid[fNEbins - 1] ? energy : fEGrid[fNEbins - 1] * 0.999;
      energy = max<double>(energy,fEGrid[0]);
      int ibin = log(energy / fEGrid[0]) * fEilDelta;
      ibin = ibin < fNEbins - 1 ? ibin : fNEbins - 2;

      double en1 = fEGrid[ibin];
      double en2 = fEGrid[ibin + 1];
      // 1.
      // this material/mixture (TMXsec object) is composed of fNElems elements
      int iel = -1; // index of elements in TEXsec** fElems  ([fNElems]) of this
      if (fNElems == 1) {
        iel = 0;
      } else { // sampling one element based on the elemntal realtive X-secs that
        // have been normalized in CTR at the Ebins!; energy interp. is
        // included now-> while loop is to avoid problems from interp.
        double xrat = (en2 - energy) / (en2 - en1);
        double xnorm = 1.;
        while (iel < 0) {
          double ran = xnorm * td->fRndm->Rndm();
          double xsum = 0;
          for (int i = 0; i < fNElems; ++i) { // simple sampling from discrete p.
            xsum += xrat * fRelXS[ibin * fNElems + i] + (1 - xrat) * fRelXS[(ibin + 1) * fNElems + i];
            if (ran <= xsum) {
              iel = i;
              break;
            }
          }
          xnorm = xsum;
        }
      }
      // at this point the index of the element is sampled:= iel
      // the corresponding TEXsec* is fElems[iel]

      // sample the reaction by using the TEXsec* that corresponds to the iel-th
      // element i.e. fElems[iel] for the current particle (with particle index of
      // ipart) at the current energy; will retrun with the reaction index;
      tracksin.fProcessV[t] = fElems[iel]->SampleReac(ipart, energy, rndArray[1]);
      tracksin.fEindexV[t] = fElems[iel]->Index(); // index of the selected element in TTabPhysMrg::fElemXsec[]
    }
  } else if (fDecayTable->HasDecay(ipart)) {
    // only decay can happen because ipart>nParticleWithReaction
    tracksin.fProcessV[t] = 3;  // decay
    tracksin.fEindexV[t] = -1;  // on nothing
  } else {                      // nothing happens
    tracksin.fProcessV[t] = -1; // nothing
    tracksin.fEindexV[t] = -1;  // on nothing
  }
}

// sample one of the elements based on #atoms/volue
//______________________________________________________________________________
int TMXsec::SampleElement(GeantTaskData *td) {
  if (fNElems > 1) {
    double randn = td->fRndm->Rndm();
    for (int itr = 0; itr < fNElems; ++itr)
      if (fRatios[itr] > randn)
        return fElems[itr]->Index(); // TTabPhysMgr index of the sampled elemet
  }

  return fElems[0]->Index();
}

// sample one of the elements based on #atoms/volue
// (same as above but with simple rng for Geant4)
//______________________________________________________________________________
int TMXsec::SampleElement() {
  if (fNElems > 1) {
    double randn = gRandom->Rndm();
    for (int itr = 0; itr < fNElems; ++itr)
      if (fRatios[itr] > randn)
        return fElems[itr]->Index(); // TTabPhysMgr index of the sampled elemet
  }

  return fElems[0]->Index();
}

//____________________________________________________________________________
int TMXsec::SelectElement(int pindex, int rindex, double energy) {
  // this material/mixture (TMXsec object) is composed of fNElems elements
  // iel is the index of elements in TEXsec** fElems  ([fNElems]) for
  // a particuclar particle type and reaction.
  // Then it returns fElems[iel]->Index()

  int iel = fNElems - 1;

  if (fNElems > 1) {
    double totalxs = 0.;
    for (int i = 0; i < fNElems; ++i) {
      // fRatios[i] is the relative #i-th atoms/volume; fRatios[] is normalized to 1.
      totalxs += fElems[i]->XS(pindex, rindex, energy) * fRatios[i];
    }

    double ranxs = totalxs * gRandom->Rndm();
    double cross = 0.;

    for (int i = 0; i < fNElems - 1; ++i) {
      // redundant, should be stored in a temporary array when calcuating totalxs
      cross += fElems[i]->XS(pindex, rindex, energy);
      if (ranxs < cross) {
        iel = i;
        break;
      }
    }
  }
  // index of the selected element
  return fElems[iel]->Index();
}

//____________________________________________________________________________
void TMXsec::Print(const char *) const {
  printf("Material %s %s with %d elements\n", GetName(), GetTitle(), fNElems);
  for (int i = 0; i < fNElems; ++i) {
    printf("%s %s\n", fElems[i]->GetName(), fElems[i]->GetTitle());
  }
}
