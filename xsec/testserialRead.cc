#include "TEXsec.h"
#include "TEFstate.h"
#include "TPDecay.h"

#include "TTabPhysMgr.h"
#include "GeantPropagator.h"

#include <iostream>
#include <fstream>

void expandPhysics(char *buf);

using std::cout;
using std::endl;

int main()
{

   char *buf=nullptr;
   int totsize;
   // read from file
   std::ifstream fin("xfphys.bin", std::ios::binary);
   fin.read(reinterpret_cast<char*>(&totsize), sizeof(totsize));
   //buf = new char[totsize];
   buf = (char*)_mm_malloc(totsize,sizeof(double));
   fin.read(reinterpret_cast<char*>(buf), totsize);
   fin.close();
   std::cout << "Total size of store " << totsize << std::endl;

   expandPhysics(buf);
  //cout<<"MAIN MAGIC "<<TEXsec::Element(10)->GetMagic()<<endl<<std::flush;
  // cout<<"MAIN XS "<< TEXsec::Element(10)->XS(10,4, 1.)<<endl<<std::flush;
   const char *fxsec = "/dev/null";
   const char *ffins = "/dev/null";
   TTabPhysMgr::Instance(fxsec, ffins );

   constexpr int nrep = 1000;
   srand(12345);
   std::ofstream fftest("xphysR.txt");
   for(auto iel=0; iel<TEXsec::NLdElems(); ++iel) {
      for(auto irep=0; irep<nrep; ++irep) {
	 // Get a random particle & reaction & energy
	 int ipart = (((double) rand())/RAND_MAX) * TPartIndex::I()->NPartReac();
//	 if (abs(TPartIndex::I()->PDG(ipart)) > 100000000) continue;
         int ireac = (((double) rand())/RAND_MAX) * FNPROC;
	 float en = (((double) rand())/RAND_MAX) * (TPartIndex::I()->Emax() - TPartIndex::I()->Emin())
	    + TPartIndex::I()->Emin();
         float xs = TEXsec::Element(iel)->XS(ipart, ireac, en);
 //        cout<<"MAIN 0 "<<iel<<" "<<irep<<endl;
 	 if(xs < 0) continue;
	 int npart=0;
	 float weight=0;
	 float kerma=0;
	 float enr=0;
	 const int *pid=0;
	 const float *mom=0;
	 int ebinindx=0;
        // cout<<" MAIN "<<ipart<<" "<<ireac<<" "<<en<<endl;
	 TEFstate::Element(iel)->SampleReac(ipart, ireac, en, npart, weight, kerma, enr, pid, mom, ebinindx);
	 if(npart <= 0) continue;
	 fftest <<  iel << ":" << TPartIndex::I()->PartName(ipart) << ":" << ireac << ":" << en
		<< ":" << xs << ":" << npart << ":" << weight << ":" << kerma << ":" << enr << ":";
	 for(auto i=0; i<npart; ++i)
	    fftest << pid[i] << ":" << mom[i*3] << ":" << mom[i*3+1] << ":" << mom[i*3+2];
	 fftest <<":" << ebinindx << std::endl;
      }
   }
   fftest.close();
   return 0;
}

void expandPhysics(char *buf) {
   std::cout << "Rebuilding TPartIndex store" << std::endl;
   TPartIndex::I()->RebuildClass(buf);
   int sizet = TPartIndex::I()->SizeOf();
   std::cout << "Number of bytes for TPartIndex " << sizet << std::endl;
   buf += sizet;
   std::cout << "Rebuilding x-sec store" << std::endl;
   TEXsec::RebuildStore(buf);
   int sizex = TEXsec::SizeOfStore();
   std::cout << "Number of bytes for x-sec " << sizex << std::endl;
   buf += sizex;
   std::cout << "Rebuilding decay store" << std::endl;
   TPDecay *dec = (TPDecay *) buf;
   dec->RebuildClass();
   TEFstate::SetDecayTable(dec);
   int sized = dec->SizeOf();
   std::cout << "Number of bytes for decay " << sized << std::endl;
   buf += sized;
   std::cout << "Rebuilding final state store" << std::endl;
   TEFstate::RebuildStore(buf);
   int sizef = TEFstate::SizeOfStore();
   std::cout << "Number of bytes for final state " << sizef << std::endl;
}
