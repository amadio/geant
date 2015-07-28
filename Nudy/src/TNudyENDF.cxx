/*
   This is the main class to read a file in ENDF format and write a file
   in R-ENDF format

   fca 2-mai-2010
*/

#include <string.h>

#include <Riostream.h>
#include <TNudyENDF.h>
#include <TNudyEndfTape.h>
#include <TNudyEndfMat.h>
#include <TNudyEndfFile.h>
#include <TNudyEndfCont.h>
#include <TNudyEndfList.h>
#include <TNudyEndfTab1.h>
#include <TNudyEndfTab2.h>
#include <TNudyEndfINTG.h>
#include <TFile.h>
#include <TMath.h>
#include <cstdlib>

const char TNudyENDF::fkElNam[119][4]={"n","H","He","Li","Be","B","C","N","O","F","Ne",
				       "Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti",
				       "V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se",
				       "Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd",
				       "Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce",
				       "Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
				       "Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb",
				       "Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu",
				       "Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg",
				       "Bh","Hs","Mt","Ds","Rg","Cn","Uut","Uuq","Uup","Uuh","Uus","Uuo"};
   
const char TNudyENDF::fkElIso[4][2]={"","m","n","o"};

ClassImp(TNudyENDF)

//_______________________________________________________________________________
TNudyENDF::TNudyENDF() :
  fLogLev(0),
  fENDF(),
  fRENDF(NULL),
  fTape(NULL),
  fMat(NULL)
{
  //
  // Standard constructor
  //
  fLine[0]='\0';
}

//_______________________________________________________________________________
TNudyENDF::TNudyENDF(const Char_t *nFileENDF, const Char_t *nFileRENDF, Option_t *opt, UChar_t loglev) :
  fLogLev(loglev),
  fENDF(),
  fRENDF(NULL),
  fTape(NULL),
  fMat(NULL)
{

  fLine[0]='\0';
  
  // Open input stream 
  fENDF.open(nFileENDF);
  if(!fENDF.is_open())
    Fatal("ctor","Could not open input file %s",nFileENDF);

  // Open output RENDF file
  fRENDF=TFile::Open(nFileRENDF,opt);
  if(!fRENDF) 
    Fatal("ctor","Could not open output file %s",nFileRENDF);
  
  // Read to Tape Identifier
  fENDF.getline(fLine,LINLEN);
  if(loglev > 3)
    std::cout << "Opened ENDF file:" << std::endl << fLine << std::endl;

  fTape = new TNudyEndfTape(fLine,fLogLev);

}

//_______________________________________________________________________________
void TNudyENDF::Process()
{
  //
  // Process a tape
  //
  Double_t c[2];
  Int_t    nl[4];
  Int_t    mtf[4];

  Int_t &curMAT = mtf[0];
  Int_t oldMAT = 0;

  while(!fENDF.eof()) {
    fENDF.getline(fLine,LINLEN);
    if(fLogLev>10) std::cout << fLine << std::endl;

    // See what we have
    GetMTF(mtf);

    if(curMAT==-1) {
      // End of Tape, finish processing
      CheckTEND();
      break;

    } else if(curMAT!=oldMAT) {
      // Normal situation, we should have read the mat till the end
      oldMAT=curMAT;
      if(fLogLev>3) std::cout << "Material(MAT) : " << curMAT << std::endl;
      // Create new material section
      GetCONT(c,nl,mtf);      
      fMat=new TNudyEndfMat(curMAT,TMath::Nint(c[0]),c[1],nl[0],(nl[1]==1),nl[2],nl[3]);
      Process(fMat);
      // Add material section to the tape list
      fTape->AddMat(fMat);
    } else  {
      // Something went wrong, we should have processed the material to the end
      Error("Process()","Did not process Material MAT %d to the end or there are mutiple evaluations of the same material\n",oldMAT);
      oldMAT=curMAT;
      if(fLogLev>3) std::cout << "Material(MAT) :  " << curMAT << std::endl;
      // Create new material section
      GetCONT(c,nl,mtf);
      fMat=new TNudyEndfMat(curMAT,TMath::Nint(c[0]),c[1],nl[0],(nl[1]==1),nl[2],nl[3]);
      Process(fMat);
      // Add material section to the tape list
      fTape->AddMat(fMat);
    }
  }

  // Write the tape to disk
  fTape->Print();
  fTape->Write();
  fRENDF->Close();
}

//_______________________________________________________________________________
void TNudyENDF::Process(TNudyEndfMat *mat)
{
  //
  // Here we read the whole material
  //
  Double_t c[2];
  Int_t    nl[4];
  Int_t    mtf[4];

  Int_t &curMAT = mtf[0];
  Int_t &curMF  = mtf[1];
  Int_t oldMF = 0;

  Char_t name[12] = {' '};
  Int_t mZA;

  // Continue reading the header 
  // Section 1 451 is put directly in the material structure
  // READ(LIB,10)C1,C2,L1,L2,N1,N2,MAT,MF,MT,NS
  GetMTF(mtf);
  if(fLogLev>3) Info("Process","File(MF) : %d",curMF);
  
  // Material record line 2
  fENDF.getline(fLine,LINLEN);
  GetCONT(c,nl,mtf);

  mat->SetELIS(c[0]);
  mat->SetSTA(c[1]);
  mat->SetLIS(nl[0]);           // Because nl's are INT and Function expects char
  mat->SetLISO(nl[1]);

  // Now we can build the name
  mZA=mat->GetZA();

  snprintf(name,13,"%d-%s-%03d%s",mZA/1000,fkElNam[mZA/1000],mZA%1000,fkElIso[mat->GetLISO()]);
  
  mat->SetName(name);
  if(fLogLev>4) Info("Process","Material Name : %s",name);
  // End of build of material name
  if(nl[2])
    Error("Process","Format error, N1 = %d in the second record of %s\n",nl[2],mat->GetName());
  mat->SetNFOR(nl[3]);   

  // Material record line 3
  if(fLogLev>4) Info("Process","Reading Header Line 3");
  fENDF.getline(fLine,LINLEN);
  GetCONT(c,nl,mtf);

  mat->SetAWI(c[0]);
  mat->SetEMAX(c[1]);
  mat->SetLREL(nl[0]);
  if(nl[1])
    Error("Process","Format error, L2 = %d in the third record of %s\n",nl[1],mat->GetName());
  mat->SetNSUB(nl[2]);
  mat->SetNVER(nl[3]);

  // Material record line 4
  if(fLogLev>4) Info("Process","Reading Header Line 4");
  fENDF.getline(fLine,LINLEN);
  GetCONT(c,nl,mtf);
  mat->SetTEMP(c[0]);
  if(c[1])
    Error("Process","Format error, C2 = %f in the fourth record of %s\n",c[1],mat->GetName());
  mat->SetLDRV(nl[0]);
  if(nl[1])
    Error("Process","Format error, L2 = %d in the fourth record of %s\n",nl[1],mat->GetName());
  mat->SetNWD(nl[2]);
  mat->SetNXC(nl[3]);
  
  // Material record line 5
  if(fLogLev>4) Info("Process","Reading Header Line 5");
  fENDF.getline(fLine,LINLEN);
  mat->SetZSYMAM(fLine);
  mat->SetALAB(&fLine[11]);
  mat->SetEDATE(&fLine[22]);
  mat->SetAUTH(&fLine[33]);

  // Material record line 6
  if(fLogLev>4) Info("Process","Reading Header Line 6");
  fENDF.getline(fLine,LINLEN);
  mat->SetREF(&fLine[1]);
  mat->SetDDATE(&fLine[22]);
  mat->SetRDATE(&fLine[33]);
  mat->SetENDATE(atoi(&fLine[55]));

  // Getting the DESC records
  fENDF.getline(fLine,LINLEN);
  mat->SetHSUB(fLine,0);
  if(fLogLev>4) Info("Process","Reading Description Record 1");
  fENDF.getline(fLine,LINLEN);
  mat->SetHSUB(fLine,1);
  if(fLogLev>4) Info("Process","Reading Description Record 2");
  fENDF.getline(fLine,LINLEN);
  mat->SetHSUB(fLine,2);
  if(fLogLev>4) Info("Process","Reading Description Record 3");
  // Get the rest of the description recrods
  if(fLogLev>6) Info("Process","No of Description TEXT records : %d",mat->GetNWD());
  int descLen=67;
  for(Int_t i=0; i<mat->GetNWD()-5; i++) {
    Char_t desc[descLen];
    fENDF.getline(fLine,LINLEN);
    strncpy(desc,fLine,descLen-1);
    desc[descLen-1]='\0';
    mat->SetDesc(desc,i);
  }

  // Get the dictionary
  if(fLogLev>6) Info("Process","Length of the Dictionary : %d",mat->GetNXC());

  for(Int_t i=0; i<mat->GetNXC(); ++i) {
    fENDF.getline(fLine,LINLEN);
    GetCONT(c,nl,mtf);
    mat->SetMFn(nl[0],i);
    mat->SetMTn(nl[1],i);
    mat->SetNCn(nl[2],i);
    mat->SetMODn(nl[3],i);
  }

  
  // Check that the next line is the SEND of 1 451
	
  GetSEND(mtf);
  //  oldMF = -1;
  // Start normal processing of the rest of material
  while(!fENDF.eof()) {
    fENDF.getline(fLine,LINLEN);
    if(fLogLev>10) std::cout << fLine << std::endl;

    // See what we have
    GetMTF(mtf);
    if(!curMAT) {
      // End of Mat, return from the method
      CheckMEND(mtf);
      break;

    } else if(curMF!=oldMF) {
      // Normal situation, we should have read the file till the end
      // Note that MF=1 is already partially read, so we process the other sections
      oldMF=curMF;
      if(fLogLev>3 && curMF!=1) Info("Process","File(MF) : %d   Section(MT) : %d",curMF,mtf[2]);
      // Create new file
      TNudyEndfFile *file=new TNudyEndfFile(mat->GetMAT(),curMF);
      Process(file);
      // Add file section to the material list
      mat->Add(file);

    } else 
      if(curMF==0 && oldMF==0) {
	// If MF 1 has only the 451 section check that it is a FEND and go on
	TNudyEndfFile *file1=new TNudyEndfFile(mat->GetMAT(),1);
	mat->Add(file1);
	CheckFEND(mtf);
      }
      else{
	// Something went wrong, we should have processed the material to the end
	Error("Process(TNudyEndfMat*)","Did not process File MF %d of Material MAT %d to the end or multiple files of the same type exist\n",curMF,curMAT);
	TNudyEndfFile *file1=new TNudyEndfFile(mat->GetMAT(),1);
	mat->Add(file1);
	std::cout << "Into checkfend 2" << std::endl;
	CheckFEND(mtf);
      }

  }
}

//_______________________________________________________________________________
void TNudyENDF::Process(TNudyEndfFile *file)
{
  Double_t c[2];
  Int_t nl[4];
  Int_t mtf[4];
  Int_t &curMAT = mtf[0];
  Int_t &curMF  = mtf[1];
  Int_t &curMT  = mtf[2];
  Int_t oldMT   = 0;

  while(!fENDF.eof()) {
    if(fLogLev>10) std::cout << fLine << std::endl;

    // See what we have
    // Careful to the inverted logics... we are already at the beg of the new file here
    // so we cannot read the next line, otherwise we skip the first line of the new file
    GetMTF(mtf);
    
    if(!curMF) {
      // End of File, return from the method
	std::cout << "Into checkfend 3" << std::endl;
      CheckFEND(mtf);
      break;
    } else if(curMT!=oldMT) {
      // Normal situation, we should have read the section till the end
      oldMT=curMT;
      if(fLogLev>3) Info("Process","New Section(MT) : %d",curMT);
      // Create new section
      GetCONT(c,nl,mtf);      
      TNudyEndfSec *sec = new TNudyEndfSec(curMAT,curMF,curMT,c[0],c[1],nl[0],nl[1],nl[2],nl[3]);          
      Process(sec);
      // Add setion to file list
      file->Add(sec);     
    } else 
      Error("Process(TNudyEndfFile*)","Did not process Section MT %d of File MF %d to the end\n",oldMT,file->GetMF());
    fENDF.getline(fLine,LINLEN);
 }
}

//_______________________________________________________________________________
void TNudyENDF::Process(TNudyEndfSec *sec)
{
  Double_t c[2] = {0,0};
  Int_t nl[4] = {0,0,0,0};
  Int_t mtf[4];
  Int_t &curMF  = mtf[1];
  Int_t &curMT  = mtf[2];

  // We are at the beginning of a section, we get the head
  GetMTF(mtf);

  switch(curMF) {
  case 1:    // ------------------- File 1
    ProcessF1(sec);
    break;
  case 2:    // ------------------- File 2
    ProcessF2(sec);
    break;
  case 3:    // ------------------- File 3
    ProcessF3(sec);
    break;
  case 4:    // ------------------- File 4
    ProcessF4(sec);
    break;
  case 5:    // ------------------- File 5
    ProcessF5(sec);
    break;
  case 6:    // ------------------- File 6
    ProcessF6(sec);
    break;
  case 7:    // ------------------- File 7
    ProcessF7(sec);
    break;
  case 8:    // ------------------- File 8
    ProcessF8(sec);
    break;
  case 9:    // ------------------- File 9
    ProcessF9(sec);
    break;
  case 10:    // ------------------- File 10
    ProcessF10(sec);
    break;
  case 12:    // ------------------- File 12
    ProcessF12(sec);
    break;
  case 13:    // ------------------- File 13
    ProcessF13(sec);
    break;
  case 14:    // ------------------- File 14
    ProcessF14(sec);
    break;
  case 15:    // ------------------- File 15
    ProcessF15(sec);
    break;
  case 23:    // ------------------- File 23
    ProcessF23(sec);
    break;
  case 26:    // ------------------- File 26
    ProcessF26(sec);
    break;
  case 27:    // ------------------- File 27
    ProcessF27(sec);
    break;
  case 28:    // ------------------- File 28
    ProcessF28(sec);
    break;
  case 31:    // ------------------- File 31
    ProcessF31(sec);
    break;
  case 32:    // ------------------- File 32
    ProcessF32(sec);
    break;
  case 33:    // ------------------- File 33
    ProcessF33(sec);
    break;
  case 34:    // ------------------- File 34
    ProcessF34(sec);
    break;
  case 35:    // ------------------- File 35
    ProcessF35(sec);
    break;
  case 40:    // ------------------- File 40
    ProcessF40(sec);
    break;
  default:    // ------------------- Untreated section of any other file
    // Just create a CONT record section
    TNudyEndfCont *secCont = new TNudyEndfCont(c[0],c[1],
					       nl[0],nl[1],nl[2],nl[3]);
    sec->Add(secCont);
    // Read the section
    if(curMT != 451)
    Error("Process(mat)","Untreated File MF %d (Section MT %d)\n",curMF,curMT);
    ToEndSec();
    break;
  }
}

//_______________________________________________________________________________
void TNudyENDF::ProcessF1(TNudyEndfSec *sec)
{
  Double_t c[2];
  Int_t nl[4];
  Int_t mtf[4];
  Int_t &curMAT = mtf[0];
  Int_t &curMF  = mtf[1];
  Int_t &curMT  = mtf[2];

  GetCONT(c,nl,mtf);
  if(curMF!=1) 
    Fatal("ProcessF1(TNudyEndfSec*)","File %d should not be processed here",curMF);

  switch(curMT) {
  case 451:
    break;
  case 452:    // ------------------- Section 452 -- Fission Yield
    switch(nl[1]) {
    case 1:    // ------------------- Section 452 -- Fission Yield LNU = 1
      {	
	TNudyEndfList *secList = new TNudyEndfList();
	// Read the section
	Process(secList);
	sec->Add(secList);
	// Get the Section END record and check it
	std::cout << "into getsend 33" << std::endl;
	GetSEND(mtf);
      }
      break;
    case 2:    // ------------------- Section 452 -- Fission Yield LNU = 2
      {
	TNudyEndfTab1 *secTab1 = new TNudyEndfTab1();
	// Read the section
	Process(secTab1);
	sec->Add(secTab1);
	// Get the Section END record and check it
	std::cout << "into getsend 1" << std::endl;
	GetSEND(mtf);
      }
      break;
    default:
      Fatal("ProcessF1(TNudyEndfSec*)","Unknown LNU value %d for Material MT %d, File MF %d and Section MT %d",
	    nl[1], curMAT, curMT, curMF);
      break;
    }
    break;
  case 455:    // ------------------- Section 455 -- Delayed neutron data
    switch(nl[1]) {
    case 1:    // ------------------- Section 455 -- LNU=1 Polynomial expansion
    case 2:    // ------------------- Section 455 -- LNU=2 Tabulated data
      switch(nl[0]) {
      case 0:    // ------------------- Section 455 -- LNU=1&2 LDG=0 energy independent
	{
	  TNudyEndfList *secListD = new TNudyEndfList();
	  // Read the section
	  Process(secListD);
	  sec->Add(secListD);
	  TNudyEndfTab1 *secTab1 = new TNudyEndfTab1();
	  // Read the section
	  Process(secTab1);
	  sec->Add(secTab1);
	  // Get the Section END record and check it
	std::cout << "into getsend 2" << std::endl;
	  GetSEND(mtf);
	}
	break;
      case 1:    // ------------------- Section 455 -- LNU=1&2 LDG=1 energy dependent
	{
	  TNudyEndfTab2 *secTab2 = new TNudyEndfTab2();
	  // Read the section
	  Process(secTab2);
	  sec->Add(secTab2);
	  // Read all sections
	  for(Int_t i=0; i<secTab2->GetNZ(); ++i) {
	    TNudyEndfList *secList = new TNudyEndfList();
	    // Read the section
	    Process(secList);
	    sec->Add(secList);
	  }
	  TNudyEndfTab1 *secTab1 = new TNudyEndfTab1();
	  // Read the section
	  Process(secTab1);
	  sec->Add(secTab1);
	  // Get the Section END record and check it
	std::cout << "into getsend 3" << std::endl;
	  GetSEND(mtf);
	}
	break;
      default:
	Fatal("ProcessF1(TNudyEndfSec*)","Unknown LDG value %d for LNU %d for Material MT %d, File MF %d and Section MT %d",
	      nl[0],nl[1], curMAT, curMT, curMF);
	break;
      }
      break;
    default:
      Fatal("ProcessF1(TNudyEndfSec*)","Unknown LNU value %d for Material MT %d, File MF %d and Section MT %d",
	    nl[1], curMAT, curMT, curMF);
      break;
    }
    break;
  case 456:    // ------------------- Section 456 -- Number of Prompt Neutrons per Fission
    switch(nl[1]) {
    case 1:    // ------------------- Section 456 -- LNU=1 polynomial representation
      {
	TNudyEndfList *secList = new TNudyEndfList();
	// Read the section
	Process(secList);
	sec->Add(secList);
	// Get the Section END record and check it
	std::cout << "into getsend 4" << std::endl;
	GetSEND(mtf);
      }
      break;
    case 2:    // ------------------- Section 456 -- LNU=2 Tabulated representation
      {
	TNudyEndfTab1 *secTab1 = new TNudyEndfTab1();
	// Read the section
	Process(secTab1);
	sec->Add(secTab1);
	// Get the Section END record and check it
	std::cout << "into getsend 5" << std::endl;
	GetSEND(mtf);
      }
      break;
    default:
      Fatal("ProcessF1(TNudyEndfSec*)","Unknown LNU value %d for Material MT %d, File MF %d and Section MT %d",
	    nl[1], curMAT, curMT, curMF);
      break;
    }
    break;
  case 458:    // ------------------- Section 458 -- Components of Energy Release Due to Fission
    {
      TNudyEndfList *secList = new TNudyEndfList();
      // Read the section
      Process(secList);
      sec->Add(secList);
      // Get the Section END record and check it
	std::cout << "into getsend 6" << std::endl;
      GetSEND(mtf);
    }
    break;

  case 460:
    {
      Int_t iLO = sec->GetL1();
      if(iLO == 1){
	Int_t iNG = sec->GetN1();
	for(Int_t i = 0; i < iNG; i++){
	  TNudyEndfTab1 *secTabPD = new TNudyEndfTab1();
	  Process(secTabPD);
	  sec->Add(secTabPD);
	}
      }
      else if (iLO == 2){
	TNudyEndfList *secListPD = new TNudyEndfList();
	Process(secListPD);
	sec->Add(secListPD);
      }
      else{
	Error("ProcessF1(TNudyEndfSec*)","Invalid value of LO=%d in File %d Section %d\n",iLO,curMF,curMT);
      }
	std::cout << "into getsend 7" << std::endl;
      GetSEND(mtf);
    }
    break;
  default:    // ------------------- Untreated section of file 1
    Fatal("ProcessF1(TNudyEndfSec*)","Untreated File MF %d Section MT %d\n",curMF,curMT);
    break;
  }
}

//_______________________________________________________________________________
void TNudyENDF::ProcessF2(TNudyEndfSec *sec)
{
  Double_t c[2];
  Int_t nl[4];
  Int_t mtf[4];
  Int_t &curMF  = mtf[1];
  Int_t &curMT  = mtf[2];
  Int_t  iLRU   = -1;
  Int_t  iLRF   = -1;
  Int_t  iNRO   = -1;
  Int_t iNAPS   = -1;
  Int_t  iLFW   = -1;

  GetCONT(c,nl,mtf);

  if(curMF!=2) 
    Fatal("ProcessF2(TNudyEndfSec*)","File %d should not be processed here",curMF);
  if(fMat->GetLRP()==-1) 
    Fatal("ProcessF2(TNudyEndfSec*)","File %d should not exist!",curMF);

  switch(curMT) {
  case 151:    // ------------------- The only section of file 2
    {
      
      // Loop over isotopes
      for(Int_t i=0; i<sec->GetN1(); ++i) {
	TNudyEndfCont *secContIso=new TNudyEndfCont(); // CONT (isotope) record	
        Process(secContIso);
        sec->Add(secContIso);
        iLFW = secContIso->GetL2();
	
      // Loop over energy ranges
       for(Int_t j=0; j<secContIso->GetN1(); ++j) { 
  	  //secContIso->GetN1() is NER
          TNudyEndfCont *secContRange=new TNudyEndfCont(); // CONT (range) record

          Process(secContRange);
          sec->Add(secContRange);
          iLRU = secContRange->GetL1();
          iLRF = secContRange->GetL2();
          iNRO = secContRange->GetN1();
          iNAPS = secContRange->GetN2();	
          switch(iLRU) {
            case 0: 
              {
                TNudyEndfCont *secCont=new TNudyEndfCont();
                Process(secCont);
                sec->Add(secCont); 
              }
              break;
            case 1: // Resolved resonance parameters
              {
                switch(iLRF) {
                  case 1: // Single-level Breit-Wigner
                  case 2: // Multi-level Breit-Wigner
                  case 3: // Reich-Moore
                    {
                      if (iNRO != 0) {
                        TNudyEndfTab1 *secTab1 = new TNudyEndfTab1();
                        Process(secTab1);
                        sec->Add(secTab1);
                      }
                      TNudyEndfCont *secContLval=new TNudyEndfCont(); // CONT (L-value) record
                      Process(secContLval);
                      sec->Add(secContLval);

                      for(Int_t k=0; k<secContLval->GetN1(); ++k) {
			// loop over L-values, subsection
                        // secContLval->GetN1() is NLS
                        TNudyEndfList *secList = new TNudyEndfList();
                        Process(secList);			
                        sec->Add(secList);
			sec->ls();
                      }
                    }
                    break;
                  case 4: // Adler-Adler
                    {
                      if (iNRO != 0) {
                        TNudyEndfTab1 *secTab1 = new TNudyEndfTab1();
                        Process(secTab1);
                        sec->Add(secTab1);
                      }
                      
                      TNudyEndfCont *secContLval=new TNudyEndfCont(); // CONT (L-value) record
                      Process(secContLval);
                      sec->Add(secContLval);
                      for(Int_t k=0; k<secContLval->GetN1(); ++k) {
                        // loop over L-values, subsection
                        // secContLval->GetN1() is NLS
                        TNudyEndfList *secList = new TNudyEndfList();
                        Process(secList);
                        sec->Add(secList);
                      }
                      TNudyEndfCont *secContJval=new TNudyEndfCont(); // CONT (J-value) record
                      Process(secContJval);
                      sec->Add(secContJval);
                      for(Int_t k=0; k<secContJval->GetN1(); ++k) {
                        // loop over J-values, subsection
                        // secContLval->GetN1() is NJS
                        TNudyEndfList *secList = new TNudyEndfList();
                        Process(secList);
                        sec->Add(secList);
                      }
                    }
                    break;
                  case 5: // This option is no longer available
                    Fatal("ProcessF2(TNudyEndfSec*)","Unknown values of LRF %d",iLRF);
                    break;
                  case 6: // This option is no longer available
                    Fatal("ProcessF2(TNudyEndfSec*)","Unknown values of LRF %d",iLRF);
                    break;
                  case 7: // R-Matrix Limited format
                    {
                      TNudyEndfCont *secContJval=new TNudyEndfCont(); // CONT (J-value) record
                      Process(secContJval);
                      sec->Add(secContJval);
                      
                      TNudyEndfList *secList = new TNudyEndfList(); // particle-pair descriptions
                      Process(secList);
                      sec->Add(secList);
                      
                      for(Int_t k=0; k<secContJval->GetN1(); ++k) {
                        // loop over NJS spin group
                        TNudyEndfList *secListCh = new TNudyEndfList(); // channel description
                        Process(secListCh);
                        sec->Add(secListCh);

                        // Formats for optional extensions to the RML
                        if (secContJval->GetL2() > 0) {
                          TNudyEndfList *secListRML = new TNudyEndfList(); 
                          Process(secListRML);
                          sec->Add(secListRML);
                          if (secListRML->GetN2()==1) { // LBK = 1 option
                            TNudyEndfTab1 *secTab1RBR = new TNudyEndfTab1(); // RBR(E)
                            Process(secTab1RBR);
                            sec->Add(secTab1RBR);
                            
                            TNudyEndfTab1 *secTab1RBI = new TNudyEndfTab1(); // RBI(E)
                            Process(secTab1RBI);
                            sec->Add(secTab1RBI);
                          }
                        } // I'm not sure this formats is correct or not, this should be tested.
                        
                        TNudyEndfList *secListR = new TNudyEndfList(); // resonance energy and widths
                        Process(secListR);
                        sec->Add(secListR);

                      }
                    }
                    break;
                  default:
                    Fatal("ProcessF2(TNudyEndfSec*)","Unknown values of LRF %d",iLRF);
                    break;
                }
              }
              break;
            case 2: // Unresolved resonance parameters
              {
                switch(iLRF) {
                  case 1:
                    {
                      switch(iLFW) { 
                        case 0:
                          // LFW = 0: fission widths not given
                          // LRF = 1: all parameters are energy-independent
                          { 
                            TNudyEndfCont *secContLval=new TNudyEndfCont(); // CONT (L-value) record
                            Process(secContLval);
                            sec->Add(secContLval);
                            for(Int_t k=0; k<secContLval->GetN1(); ++k) {
                              // loop over L-values, subsection
                              // secContLval->GetN1() is NLS
                              TNudyEndfList *secList = new TNudyEndfList();
                              Process(secList);
                              sec->Add(secList);
                            }
                          }
                          break;
                        case 1:
                          // LFW = 1: fission widths given
                          // LRF = 1: only fission widths are energy-dependent
                          { 	
                            TNudyEndfList *secList = new TNudyEndfList(); // This List record controls J-value records
                            Process(secList);
                            sec->Add(secList);
                        
                            for(Int_t l=0; l<secList->GetN2(); ++l) {                        
                              TNudyEndfCont *secContJval=new TNudyEndfCont(); // CONT (J-value) record
                              Process(secContJval);
                              sec->Add(secContJval);

                              for(Int_t k=0; k<secContJval->GetN1(); ++k) {
                                // loop over J-values, subsection
                                // secContLval->GetN1() is NJS
                                TNudyEndfList *secListNJS = new TNudyEndfList();
                                Process(secListNJS);
                                sec->Add(secListNJS);
                              }
                            }
                          }
                          break;
                        default:
                          Fatal("ProcessF2(TNudyEndfSec*)","Unknown values of LFW %d",iLFW);
                          break;
                      }
                      break;
                    }
                  case 2: // LRF = 2 all energy-dependent parameters
                    {
                      TNudyEndfCont *secContLval=new TNudyEndfCont(); // CONT (L-value) record
                      Process(secContLval);
                      sec->Add(secContLval);
                      for (Int_t l=0; l<secContLval->GetN1(); ++l) {   
                        TNudyEndfCont *secContJval=new TNudyEndfCont(); // CONT (J-value) record
                        Process(secContJval);
                        sec->Add(secContJval);

                        for(Int_t k=0; k<secContJval->GetN1(); ++k) {
                          // loop over J-values, subsection
                          // secContLval->GetN1() is NJS
                          TNudyEndfList *secList = new TNudyEndfList();
                          Process(secList);
                          sec->Add(secList);
                        }
                      }				  
                    }
                    break;
                  default:
                    Fatal("ProcessF2(TNudyEndfSec*)","Unknown values of LRF %d",iLRF);
                    break;
                }
	      }
              break;
            default:
              Fatal("ProcessF2(TNudyEndfSec*)","Unknown values of LRU %d",iLRU);
              break;
          }	
        }
      }
    }
    break;
  default:    // ------------------- Untreated section of file 2
    Fatal("ProcessF2(TNudyEndfSec*)","Untreated File MF %d Section MT %d\n",curMF,curMT);
    break;
  }
	std::cout << "into getsend 8" << std::endl;
  GetSEND(mtf); 
}
//_______________________________________________________________________________
void TNudyENDF::ProcessF3(TNudyEndfSec *sec)
{
  Double_t c[2];
  Int_t  nl[4];
  Int_t mtf[4];
  Int_t &curMF  = mtf[1];
  
  GetCONT(c,nl,mtf);
  if(curMF!=3) 
    Fatal("ProcessF3(TNudyEndfSec*)","File %d should not be processed here",curMF);
  
  TNudyEndfTab1 *secTab1 = new TNudyEndfTab1();
  // Read the section
  Process(secTab1);
  sec->Add(secTab1);
  // Get the Section END record and check it
	std::cout << "into getsend 9" << std::endl;
  GetSEND(mtf);
}

//_______________________________________________________________________________
void TNudyENDF::ProcessF4(TNudyEndfSec *sec)
{
  Double_t c[2];
  Int_t  nl[4];
  Int_t mtf[4];
  Int_t &curMF  = mtf[1];
  Int_t iLTT = -1;
  Int_t  iLI = -1;
    
  GetCONT(c,nl,mtf);
  if(curMF!=4) 
    Fatal("ProcessF4(TNudyEndfSec*)","File %d should not be processed here",curMF);
  if(sec->GetL1()!=0) {  // Old format, this is not available anymore
    Error("ProcessF4(TNudyEndfSec*)","LVT %d is not available, it is old format", sec->GetL1());
    ToEndSec();
  }
  else {  // Normal format
    iLTT = sec->GetL2();
    
    TNudyEndfCont *secCont=new TNudyEndfCont(); // CONT (representation) record
    Process(secCont);
    sec->Add(secCont);
    iLI = secCont->GetL1();
    
    switch(iLTT) {
      case 0: {    // LTT=0, LI=1: Purely Isotropic Angular Distributions
        if (iLI!=1)
          Fatal("ProcessF4(TNudyEndfSec*)","Unknown values of LI %d",iLI);
      }
      break;
      case 1: {    // LTT=1, LI=0: Legendre Polynomial Coefficients
        TNudyEndfTab2 *secTab2 = new TNudyEndfTab2();
  	    Process(secTab2); 
  	    sec->Add(secTab2);
        for(Int_t k=0; k<secTab2->GetNZ(); k++) {
          TNudyEndfList *secList = new TNudyEndfList();
          Process(secList);
          sec->Add(secList);
        }        
      }
      break;
      case 2: {    // LTT=2, LI=0: Tabulated Probability Distributions

        TNudyEndfTab2 *secTab2 = new TNudyEndfTab2();
  	    Process(secTab2);
  	    sec->Add(secTab2);
      
        for(Int_t k=0; k<secTab2->GetNZ();k++) {
          TNudyEndfTab1 *secTab1 = new TNudyEndfTab1();
          Process(secTab1);
          sec->Add(secTab1);
        } 
      }
      break;
      case 3: {
        TNudyEndfTab2 *secTab2Legen = new TNudyEndfTab2(); // Legendre
	      Process(secTab2Legen);
	      sec->Add(secTab2Legen);
        
        for(Int_t k=0; k<secTab2Legen->GetNZ(); k++) {
          TNudyEndfList *secList = new TNudyEndfList();
          Process(secList);
          sec->Add(secList);
        }  

        TNudyEndfTab2 *secTab2Tab = new TNudyEndfTab2();  // Tabulated
	      Process(secTab2Tab);
	      sec->Add(secTab2Tab);
        
        for(Int_t k=0; k<secTab2Tab->GetNZ(); k++) {
          TNudyEndfTab1 *secTab1 = new TNudyEndfTab1();
          Process(secTab1);
          sec->Add(secTab1);
        } 
      }
      break;
      default:
        Fatal("ProcessF4(TNudyEndfSec*)","Unknown values of LTT %d",iLTT);
        break;
    }
    GetSEND(mtf);
  }
}

//_______________________________________________________________________________
void TNudyENDF::ProcessF5(TNudyEndfSec *sec)
{
  Double_t c[2];
  Int_t  nl[4];
  Int_t mtf[4];
  Int_t &curMF  = mtf[1];
  Int_t iNK = -1;
  Int_t iLF = -1;

  GetCONT(c,nl,mtf);
  if(curMF!=5) 
    Fatal("ProcessF5(TNudyEndfSec*)","File %d should not be processed here",curMF);

  iNK = sec->GetN1();

  for(Int_t i=0; i<iNK; ++i) {
    //loop over partial energy distributions

    TNudyEndfTab1 *secTab1 = new TNudyEndfTab1(); // p(E) tab1
    Process(secTab1);
    sec->Add(secTab1);
    iLF = secTab1->GetL2();

    switch(iLF) {
      case  1:  // Arbitrary Tabulated Function, LF=1
        {
          TNudyEndfTab2 *secTab2 = new TNudyEndfTab2();  // Tabulated
	        Process(secTab2);
	        sec->Add(secTab2);
      
          for(Int_t k=0; k<secTab2->GetNZ(); ++k) {
            // secTab2->GetNZ() is number of incident energy
            TNudyEndfTab1 *secTab1NZ = new TNudyEndfTab1();
            Process(secTab1NZ);
            sec->Add(secTab1NZ);
          } 
        }
      break;
      case  5:  // General Evaporation Spectrum, LF=5
      case  7:  // Simple Maxwellian Fission Spectrum, LF=7
      case  9:  // Evaporation Spectrum, LF=9
      case 11:  // Energy-Dependent Watt Spectrum, LF=11
      case 12:  // Energy-Dependent Fission Neutron Spectrum, LF=12
        {
          TNudyEndfTab1 *secTab1A = new TNudyEndfTab1(); 
          // theta(E) for LF = 5,7,9
          // a(E) for LF = 11
          // TM(E) for LF = 12
          Process(secTab1A);
          sec->Add(secTab1A);
          
          if (iLF==5 || iLF==11) {// LF=5,11 have two Tab1 records 
            TNudyEndfTab1 *secTab1B = new TNudyEndfTab1(); 
            // g(x) for LF = 5
            // b(E) for LF = 11
            Process(secTab1B);
            sec->Add(secTab1B); 
          }
        }
      break;
      default:
        Fatal("ProcessF5(TNudyEndfSec*)","Unknown values of LF %d",iLF);
        break;
    }
  }
  GetSEND(mtf);
}

//_______________________________________________________________________________
void TNudyENDF::ProcessF6(TNudyEndfSec *sec)
{
  Double_t c[2];
  Int_t  nl[4];
  Int_t mtf[4];
  Int_t &curMF  = mtf[1];
  Int_t iNK     = -1;
  Int_t iLAW    = -1;

  GetCONT(c,nl,mtf);
  if(curMF!=6) 
    Fatal("ProcessF6(TNudyEndfSec*)","File %d should not be processed here",curMF);

  iNK = sec->GetN1();

  for(Int_t i=0; i<iNK; ++i) {
    // NK is number of subsections
    TNudyEndfTab1 *secTab1 = new TNudyEndfTab1(); 
    Process(secTab1);
    sec->Add(secTab1);
    iLAW = secTab1->GetL2();

    switch(iLAW){
      case 0: // unknown distribution
      case 3: // isotropic two-body distribution
      case 4: // recoil distribution of a two-body reaction
        break;
      case 1: // continuum energy-angle distribution
      case 2: // two-body reaction angular distribution
      case 5: // charged-particle elastic distribution
        {
        
          TNudyEndfTab2 *secTab2 = new TNudyEndfTab2();  
          // normal Tab2 interpolation parameters 
	        Process(secTab2);
	        sec->Add(secTab2);
          
          for(Int_t j=0; j<secTab2->GetN2(); ++j) {
            // secTab2->GetN2() is NE
            TNudyEndfList *secList = new TNudyEndfList();
            Process(secList);
            sec->Add(secList);
          }
        }
        break;
      case 6: {
        // n-body phase-space distribution
          TNudyEndfCont *secCont=new TNudyEndfCont(); // CONT record
          Process(secCont);
          sec->Add(secCont);
        }
        break;
      case 7: {
        // laboratory angle-energy law
          TNudyEndfTab2 *secTab2E = new TNudyEndfTab2();  
          // normal Tab2 interpolation parameters for incident energy, E
	        Process(secTab2E);
	        sec->Add(secTab2E);
          
          for(Int_t j=0; j<secTab2E->GetN2(); ++j) {
            // secTab2E->GetN2() is NE
            TNudyEndfTab2 *secTab2Mu = new TNudyEndfTab2();  
            // normal Tab2 interpolation parameters for emission cosine, Mu
	          Process(secTab2Mu);
	          sec->Add(secTab2Mu);
            
            for(Int_t k=0; k<secTab2Mu->GetN2(); ++k) {
              // secTab2Mu->GetN2() is NMU
              TNudyEndfTab1 *secTab1C = new TNudyEndfTab1(); 
              Process(secTab1C);
              sec->Add(secTab1C);
            }
          }
        }
        break;
      default:
        Fatal("ProcessF6(TNudyEndfSec*)","Unknown values of LAW %d",iLAW);
        break;
    }
  }
  GetSEND(mtf);
}

//_______________________________________________________________________________
void TNudyENDF::ProcessF7(TNudyEndfSec *sec)
{
  Double_t c[2];
  Int_t  nl[4];
  Int_t mtf[4];
  Int_t &curMF  = mtf[1];
  Int_t &curMT  = mtf[2];

  GetCONT(c,nl,mtf);
  if(curMF!=7) 
    Fatal("ProcessF7(TNudyEndfSec*)","File %d should not be processed here",curMF);

  switch(curMT) {
    case 2: // Elastic scattering (MT = 2)
      {
        switch(sec->GetL1()) { //sec->GetL1() is LTHR
          case 1: // Coherent elastic scattering
            {
              TNudyEndfTab1 *secTab1 = new TNudyEndfTab1(); 
              Process(secTab1);
              sec->Add(secTab1);
              
              for(Int_t i=0; i<secTab1->GetL1(); ++i) {
                TNudyEndfList *secList = new TNudyEndfList();
                Process(secList);
                sec->Add(secList);
              }
            }
            break;
          case 2: // Incoherent elastic scattering
            {
              TNudyEndfTab1 *secTab1 = new TNudyEndfTab1(); 
              Process(secTab1);
              sec->Add(secTab1);
            }
            break;
          default:
            Fatal("ProcessF7(TNudyEndfSec*)","Unknown values of LTHR %d",sec->GetL1());
            break; 
        }
      }
      break;
    case 4: // Inelastic scattering (MT = 4)
      {
        TNudyEndfList *secList = new TNudyEndfList();
        Process(secList);
        sec->Add(secList);
        TNudyEndfTab2 *secTab2 = new TNudyEndfTab2();  
        Process(secTab2);
	sec->Add(secTab2);
        for(Int_t i=0; i<secTab2->GetN2(); ++i) {
          // secTab2->GetN2() is NB(number of beta values)
          TNudyEndfTab1 *secTab1 = new TNudyEndfTab1(); 
          Process(secTab1);
          sec->Add(secTab1);
          for(Int_t j=0; j<secTab1->GetL1(); ++j) {
            // secTab1->GetL1() is LT
            // LT is the temperature dependence flag
            TNudyEndfList *secListLT = new TNudyEndfList();
            Process(secListLT);
            sec->Add(secListLT);
          }
        }
	TNudyEndfTab1 *secTabTeff0 = new TNudyEndfTab1();
	// table of effective remperatures (K)
	Process(secTabTeff0);
	sec->Add(secTabTeff0);
	if(secList->GetNPL() >= 7 && secList->GetLIST(6) == 0) {
	  TNudyEndfTab1 *secTabTeff1 = new TNudyEndfTab1();
	  // table of effective remperatures (K)
	  Process(secTabTeff1);
	  sec->Add(secTabTeff1);
	}
	if(secList->GetNPL() >= 13 && secList->GetLIST(12) == 0) {
	  TNudyEndfTab1 *secTabTeff2 = new TNudyEndfTab1();
	  // table of effective remperatures (K)
	  Process(secTabTeff2);
	  sec->Add(secTabTeff2);
	}
	if(secList->GetNPL() >= 19 && secList->GetLIST(18) == 0) {
	  TNudyEndfTab1 *secTabTeff3 = new TNudyEndfTab1();
	  // table of effective remperatures (K)
	  Process(secTabTeff3);
	  sec->Add(secTabTeff3);
	}
      }
      break;
    default:
      Fatal("ProcessF7(TNudyEndfSec*)","Untreated File MF %d Section MT %d\n",curMF,curMT);
      break;
  }
  GetSEND(mtf);
}

//_______________________________________________________________________________
void TNudyENDF::ProcessF8(TNudyEndfSec *sec)
{
  Double_t c[2];
  Int_t  nl[4];
  Int_t mtf[4];
  Int_t &curMF  = mtf[1];
  Int_t &curMT  = mtf[2];

  GetCONT(c,nl,mtf);
  if(curMF!=8) 
    Fatal("ProcessF8(TNudyEndfSec*)","File %d should not be processed here",curMF);
  
  switch(curMT) {
  case 454: // Independent fission product yield data
  case 459: // Cumulative fission product yield data
    {
      for(Int_t i=0; i<sec->GetL1(); i++) {
	TNudyEndfList *secList = new TNudyEndfList();
	Process(secList);
	sec->Add(secList);
      }
    }
    break;
  case 457: // Radioactive decay data
      {
        if(sec->GetN1()==0) {
          // for a radioactive nucleus NST=0
          TNudyEndfList *secListGI = new TNudyEndfList(); //General Information
          Process(secListGI);
          sec->Add(secListGI);

          TNudyEndfList *secListDI = new TNudyEndfList(); //Decay mode Information
          Process(secListDI);
          sec->Add(secListDI);

          for(Int_t i=0; i<sec->GetN2(); ++i) {
            TNudyEndfList *secListRRS = new TNudyEndfList(); //Resulting Radiation Spectra
            Process(secListRRS);
            sec->Add(secListRRS);
            //Editing

            if(secListRRS->GetL1() != 1) {

	      //Loop NER = secListRRS->GetN2() times
	      for(Int_t j = 0; j < secListRRS->GetN2(); j++){
		  TNudyEndfList *secListA = new TNudyEndfList(); 
		  Process(secListA);
		  sec->Add(secListA);
		}
            }
            if(secListRRS->GetL1() != 0) {
              TNudyEndfTab1 *secTab1 = new TNudyEndfTab1(); 
              Process(secTab1);
              sec->Add(secTab1);
	    
              if(secTab1->GetL2() != 0) {
                TNudyEndfList *secListB = new TNudyEndfList(); 
                Process(secListB);
                sec->Add(secListB);
              }
            }
          }
        }
        else if(sec->GetN1()==1) {
          // for a stable nucleus NST=1
          for(Int_t i=0; i<2; ++i) {
            TNudyEndfList *secListC = new TNudyEndfList();
            Process(secListC);
            sec->Add(secListC);
          }
        }
	      else
          Fatal("ProcessF8(TNudyEndfSec*)","Unknown values of NST %d",sec->GetN1());
      }
      break;
    default:
      {
        if(sec->GetN2() == 0) {
          // sec->GetN2() is NO: Flag denoting where the decay information is
          // to be given for an important radioactive end product
          // NO=0: complete decay chain given under this MT
          for(Int_t i=0; i<sec->GetN1(); ++i) {
            TNudyEndfList *secList = new TNudyEndfList();
            Process(secList);
            sec->Add(secList);
          }
        }
        else if(sec->GetN2() == 1) {
          // NO=1: decay chain given in MT=457 in the decay data file
          for(Int_t i=0; i<sec->GetN1(); ++i) {
            TNudyEndfCont *secCont=new TNudyEndfCont(); // CONT record
            Process(secCont);
            sec->Add(secCont);
          }
        }
        else
          Fatal("ProcessF8(TNudyEndfSec*)","Unknown values of NO %d",sec->GetN2());
      }
      break;
    }
  GetSEND(mtf);
}

//_______________________________________________________________________________
void TNudyENDF::ProcessF9(TNudyEndfSec *sec)
{
  Double_t c[2];
  Int_t  nl[4];
  Int_t mtf[4];
  Int_t &curMF  = mtf[1];

  GetCONT(c,nl,mtf);
  if(curMF!=9) 
    Fatal("ProcessF9(TNudyEndfSec*)","File %d should not be processed here",curMF);

  for(Int_t i=0; i<sec->GetN1(); ++i) {
    TNudyEndfTab1 *secTab1 = new TNudyEndfTab1(); 
    Process(secTab1);
    sec->Add(secTab1);
  }
  GetSEND(mtf);
}

//_______________________________________________________________________________
void TNudyENDF::ProcessF10(TNudyEndfSec *sec)
{
  Double_t c[2];
  Int_t  nl[4];
  Int_t mtf[4];
  Int_t &curMF  = mtf[1];

  GetCONT(c,nl,mtf);
  if(curMF!=10) 
    Fatal("ProcessF10(TNudyEndfSec*)","File %d should not be processed here",curMF);

  for(Int_t i=0; i<sec->GetN1(); ++i) {
    TNudyEndfTab1 *secTab1 = new TNudyEndfTab1(); 
    Process(secTab1);
    sec->Add(secTab1);
  }
  GetSEND(mtf);
}

//_______________________________________________________________________________
void TNudyENDF::ProcessF12(TNudyEndfSec *sec)
{
  Double_t c[2];
  Int_t  nl[4];
  Int_t mtf[4];
  Int_t &curMF  = mtf[1];
  Int_t iLO = -1;

  GetCONT(c,nl,mtf);
  if(curMF!=12) 
    Fatal("ProcessF12(TNudyEndfSec*)","File %d should not be processed here",curMF);

  iLO = sec->GetL1();
  switch(iLO) {
    case 1: // option 1: Multiplicities
      {
        TNudyEndfTab1 *secTab1 = new TNudyEndfTab1(); 
        Process(secTab1);
        sec->Add(secTab1);
        if(sec->GetN1()>1) for(Int_t i=0; i<sec->GetN1(); ++i) {
          TNudyEndfTab1 *secTab1A = new TNudyEndfTab1(); 
          Process(secTab1A);
          sec->Add(secTab1A);
        }
      }
      break;
    case 2: // option 2: Transition Probability Arrays
      {
        TNudyEndfList *secListB = new TNudyEndfList();
        Process(secListB);
        sec->Add(secListB);
      }
      break;
    default:
      Fatal("ProcessF12(TNudyEndfSec*)","Unknown values of LO %d",iLO);
      break;
  }
  GetSEND(mtf);
}

//_______________________________________________________________________________
void TNudyENDF::ProcessF13(TNudyEndfSec *sec)
{
  Double_t c[2];
  Int_t  nl[4];
  Int_t mtf[4];
  Int_t &curMF  = mtf[1];

  GetCONT(c,nl,mtf);
  if(curMF!=13) 
    Fatal("ProcessF13(TNudyEndfSec*)","File %d should not be processed here",curMF);

  if(sec->GetN1()!=1) {
    TNudyEndfTab1 *secTab1 = new TNudyEndfTab1(); 
    Process(secTab1);
    sec->Add(secTab1);
  }
  for(Int_t i=0; i<sec->GetN1(); ++i) {
    TNudyEndfTab1 *secTab1 = new TNudyEndfTab1(); 
    Process(secTab1);
    sec->Add(secTab1);
  }
  GetSEND(mtf);
}

//_______________________________________________________________________________
void TNudyENDF::ProcessF14(TNudyEndfSec *sec)
{
  Double_t c[2];
  Int_t  nl[4];
  Int_t mtf[4];
  Int_t &curMF  = mtf[1];
  Int_t iLI  = -1;
  Int_t iLTT = -1;
  Int_t iNK  = -1;
  Int_t iNI  = -1;

  GetCONT(c,nl,mtf);
  if(curMF!=14) 
    Fatal("ProcessF14(TNudyEndfSec*)","File %d should not be processed here",curMF);

  iLI = sec->GetL1();
  iLTT= sec->GetL2();
  iNK = sec->GetN1();
  iNI = sec->GetN2();

  switch(iLI) {
    case 1: // Isotropic Distribution
      break;
    case 0: // Anisotropic Distribution 
      {
        switch(iLTT) {
          case 1: // Legendre Coefficient Representation
            {
              if(iNI>0) for(Int_t i=0; i<iNI; ++i) {
                TNudyEndfCont *secCont=new TNudyEndfCont(); 
                Process(secCont);
                sec->Add(secCont);
              }
              for(Int_t j=0; j<(iNK-iNI); ++j) {
                TNudyEndfTab2 *secTab2 = new TNudyEndfTab2();  
                Process(secTab2);
	              sec->Add(secTab2);

                for(Int_t k=0; k<secTab2->GetN2(); ++k) {
                  TNudyEndfList *secList = new TNudyEndfList();
                  Process(secList);
                  sec->Add(secList);
                }
              }
            }
            break;
          case 2: // Tabulated Angular Distributions
            {
              if(iNI>0) for(Int_t i=0; i<iNI; ++i) {
                TNudyEndfCont *secCont=new TNudyEndfCont(); 
                Process(secCont);
                sec->Add(secCont);
              }
              for(Int_t j=0; j<(iNK-iNI); ++j) {
                TNudyEndfTab2 *secTab2 = new TNudyEndfTab2();  
                Process(secTab2);
	              sec->Add(secTab2);

                for(Int_t k=0; k<secTab2->GetN2(); ++k) {
                  TNudyEndfTab1 *secTab1 = new TNudyEndfTab1(); 
                  Process(secTab1);
                  sec->Add(secTab1);
                }
              }
            }
            break;
          default:
            Fatal("ProcessF14(TNudyEndfSec*)","Unknown values of LTT %d",iLTT);
            break;
        }
      }
      break;
    default:
      Fatal("ProcessF14(TNudyEndfSec*)","Unknown values of LI %d",iLI);
      break;
  }
  GetSEND(mtf);
}

//_______________________________________________________________________________
void TNudyENDF::ProcessF15(TNudyEndfSec *sec)
{
  Double_t c[2];
  Int_t  nl[4];
  Int_t mtf[4];
  Int_t &curMF  = mtf[1];
  Int_t iNC = -1;

  GetCONT(c,nl,mtf);
  if(curMF!=15) 
    Fatal("ProcessF15(TNudyEndfSec*)","File %d should not be processed here",curMF);

  iNC = sec->GetN1();

  for(Int_t i=0; i<iNC; ++i) {
    TNudyEndfTab1 *secTab1 = new TNudyEndfTab1(); 
    Process(secTab1);
    sec->Add(secTab1);

    if(secTab1->GetL2()!=1)
      Fatal("ProcessF15(TNudyEndfSec*)","Unknown values of LF %d",secTab1->GetL2());

    TNudyEndfTab2 *secTab2 = new TNudyEndfTab2();  
    Process(secTab2);
	  sec->Add(secTab2);

    for(Int_t j=0; j<secTab2->GetN2(); ++j) {
      TNudyEndfTab1 *secTab1A = new TNudyEndfTab1(); 
      Process(secTab1A);
      sec->Add(secTab1A);
    }
  }
  GetSEND(mtf);
}

//_______________________________________________________________________________
void TNudyENDF::ProcessF23(TNudyEndfSec *sec)
{
  Double_t c[2];
  Int_t  nl[4];
  Int_t mtf[4];
  Int_t &curMF  = mtf[1];

  GetCONT(c,nl,mtf);
  if(curMF!=23) 
    Fatal("ProcessF23(TNudyEndfSec*)","File %d should not be processed here",curMF);

  TNudyEndfTab1 *secTab1 = new TNudyEndfTab1(); 
  Process(secTab1);
  sec->Add(secTab1);
  GetSEND(mtf);
}

//_______________________________________________________________________________
void TNudyENDF::ProcessF26(TNudyEndfSec *sec)
{
  Double_t c[2];
  Int_t  nl[4];
  Int_t mtf[4];
  Int_t &curMF  = mtf[1];
  Int_t iNK     = -1;
  Int_t iLAW    = -1;

  GetCONT(c,nl,mtf);
  if(curMF!=26) 
    Fatal("ProcessF26(TNudyEndfSec*)","File %d should not be processed here",curMF);

  iNK = sec->GetN1();

  for(Int_t i=0; i<iNK; ++i) {
    // NK is number of subsections
    TNudyEndfTab1 *secTab1 = new TNudyEndfTab1(); 
    Process(secTab1);
    sec->Add(secTab1);
    iLAW = secTab1->GetL2();

    switch(iLAW){
      case 1: // continuum distribution
      case 2: // two-body angular distribution
        {        
          TNudyEndfTab2 *secTab2 = new TNudyEndfTab2();  
          // normal Tab2 interpolation parameters 
	        Process(secTab2);
	        sec->Add(secTab2);
          
          for(Int_t j=0; j<secTab2->GetN2(); ++j) {
            // secTab2->GetN2() is NE
            TNudyEndfList *secList = new TNudyEndfList();
            Process(secList);
            sec->Add(secList);
          }
        }
        break;
      case 8: // Energy transfer for excitation
        {
          TNudyEndfTab1 *secTab1A = new TNudyEndfTab1(); 
          Process(secTab1A);
          sec->Add(secTab1A);
        }
        break;
      default:
        Fatal("ProcessF26(TNudyEndfSec*)","Unknown values of LAW %d",iLAW);
        break;
    }
  }
  GetSEND(mtf);
}

//_______________________________________________________________________________
void TNudyENDF::ProcessF27(TNudyEndfSec *sec)
{
  Double_t c[2];
  Int_t  nl[4];
  Int_t mtf[4];
  Int_t &curMF  = mtf[1];

  GetCONT(c,nl,mtf);
  if(curMF!=27) 
    Fatal("ProcessF27(TNudyEndfSec*)","File %d should not be processed here",curMF);

  TNudyEndfTab1 *secTab1 = new TNudyEndfTab1(); 
  Process(secTab1);
  sec->Add(secTab1);
  GetSEND(mtf);
}

//_______________________________________________________________________________
void TNudyENDF::ProcessF28(TNudyEndfSec *sec)
{
  Double_t c[2];
  Int_t  nl[4];
  Int_t mtf[4];
  Int_t &curMF  = mtf[1];
  
  GetCONT(c,nl,mtf);
  if(curMF!=28) 
    Fatal("ProcessF28(TNudyEndfSec*)","File %d should not be processed here",curMF);
  
  for(Int_t i=0; i<sec->GetN1(); ++i) {
    TNudyEndfList *secList = new TNudyEndfList();
    Process(secList);
    sec->Add(secList);
  }
  GetSEND(mtf);
}

//_______________________________________________________________________________
void TNudyENDF::ProcessF30(TNudyEndfSec *sec)
{
  Double_t c[2];
  Int_t  nl[4];
  Int_t mtf[4];
  Int_t &curMF  = mtf[1];
  Int_t MT = -1;
  Int_t NP = -1;

  GetCONT(c,nl,mtf);
  if(curMF!=30) 
    Fatal("ProcessF30(TNudyEndfSec*)","File %d should not be processed here",curMF);

  //Read
  MT = sec->GetC1();
  NP = sec->GetN2();
  if(MT==1)
  {
	  TNudyEndfCont *secCont = new TNudyEndfCont();
	  Process(secCont);
	  sec->Add(secCont);
	  Int_t NDIR = secCont->GetN1();
	  Int_t NCTAB = secCont->GetN2();
	  for(Int_t j=0;j<NDIR;j++)
	  {
		  TNudyEndfCont *subsecCont = new TNudyEndfCont();
		  Process(subsecCont);
		  sec->Add(subsecCont);
	  }
	  for(Int_t k=0;k<NCTAB;k++)
	  {
		  TNudyEndfCont *subsecCont = new TNudyEndfCont();
		  Process(subsecCont);
		  sec->Add(subsecCont);
	  }
  }
  else if(MT==2)
  {
	for(Int_t j=0;j<NP;j++)
	{
		TNudyEndfList *subsecList = new TNudyEndfList();
		Process(subsecList);
		sec->Add(subsecList);
	}
  }
  else if(MT>=11 && MT<=999)
  {
	  Int_t NL = NP;
	  for(Int_t j=0;j<NL;j++)
	  {
		  TNudyEndfCont *subsecCont = new TNudyEndfCont();
		  Process(subsecCont);
		  sec->Add(subsecCont);
	  }
  }
  GetSEND(mtf);
}

//_______________________________________________________________________________
void TNudyENDF::ProcessF31(TNudyEndfSec *sec)
{
  Double_t c[2];
  Int_t  nl[4];
  Int_t mtf[4];
  Int_t &curMF  = mtf[1];
  //  Int_t &curMT  = mtf[2];
  Int_t &iNL = nl[3];

  GetCONT(c,nl,mtf);
  if(curMF!=31) 
    Fatal("ProcessF31(TNudyEndfSec*)","File %d should not be processed here",curMF);

  for(Int_t i=0; i<iNL; ++i) {
    TNudyEndfCont *secCont=new TNudyEndfCont(); 
    Process(secCont);
    sec->Add(secCont);
    for(Int_t j=0; j<secCont->GetN1(); ++j) {
      // NC-Type records
      TNudyEndfCont *secContNC=new TNudyEndfCont(); 
      Process(secContNC);
      sec->Add(secContNC);

      TNudyEndfList *secList = new TNudyEndfList();
      Process(secList);
      sec->Add(secList);
    }
    for(Int_t k=0; k<secCont->GetN2(); ++k) {
      // NI-Type records
      TNudyEndfList *secList = new TNudyEndfList();
      Process(secList);
      sec->Add(secList);
    }
  }
  GetSEND(mtf);
}

//_______________________________________________________________________________
void TNudyENDF::ProcessF32(TNudyEndfSec *sec)
{
  Double_t c[2];
  Int_t  nl[4];
  Int_t mtf[4];
  Int_t &curMF  = mtf[1];
  Int_t iNIS = -1;
  Int_t iLFW = -1;
  Int_t iNER = -1;
  Int_t iLRU = -1;
  Int_t iLRF = -1;
  Int_t iNRO = -1;
  Int_t iNAPS= -1;
  Int_t iLCOMP = -1;
  
  GetCONT(c,nl,mtf);
  if(curMF!=32) 
    Fatal("ProcessF32(TNudyEndfSec*)","File %d should not be processed here",curMF);

  iNIS = sec->GetN1();
  
  for(Int_t i=0; i<iNIS; ++i) {
    TNudyEndfCont *secContIso=new TNudyEndfCont(); // CONT (isotope) record
    Process(secContIso);
    sec->Add(secContIso);
    iLFW = secContIso->GetL2();
    iNER = secContIso->GetN1();
    for(Int_t j=0; j<iNER; ++j) {
      TNudyEndfCont *secContRange=new TNudyEndfCont(); // CONT (range) record
      Process(secContRange);
      sec->Add(secContRange);
      iLRU = secContRange->GetL1();
      iLRF = secContRange->GetL2();
      iNRO = secContRange->GetN1();
      iNAPS = secContRange->GetN2();
      switch(iLRU) {
      case 1: { // LRU == 1
	if(iNRO!=0){
	  // NI-Type records
	  TNudyEndfCont *secCont=new TNudyEndfCont(); 
	  Process(secCont);
	  sec->Add(secCont);
	  for(Int_t k=0; k<sec->GetN2(); ++k) {
	    TNudyEndfList *secList = new TNudyEndfList();
	    Process(secList);
	    sec->Add(secList);      
	  }
	}
	TNudyEndfCont *secContLCOMP=new TNudyEndfCont(); 
	Process(secContLCOMP);
	sec->Add(secContLCOMP);
	iLCOMP = secContLCOMP->GetL2();
	switch(iLCOMP) {
	case 0: { // LCOMP == 0
	  // Compatible resolved resonance format
	  for(Int_t k=0; k<secContLCOMP->GetN1(); ++k) {
	    TNudyEndfList *secList = new TNudyEndfList();
	    Process(secList);
	    sec->Add(secList); 
	      }
	}
	  break;
	case 1: { // LCOMP == 1
	  // General resolved resonance format
	  switch(iLRF) {
	  case 1: // Single-Level Breit-Wigner
	  case 2: // Multi-Level Breit-Wigner
	  case 4: // Adler-Adler
	    {
	      TNudyEndfCont *secCont=new TNudyEndfCont(); 
	      Process(secCont);
	      sec->Add(secCont);
	      for(Int_t l=0; l<secCont->GetN1(); ++l) {
		TNudyEndfList *secList = new TNudyEndfList();
		Process(secList);
		sec->Add(secList);
	      }
	      for(Int_t m=0; m<secCont->GetN2(); ++m) {
		TNudyEndfList *secList = new TNudyEndfList();
		Process(secList);
                sec->Add(secList);
	      }
                      }
	    break;
	  case 3: // Reich-Moore
	    {
	      TNudyEndfCont *secContNL=new TNudyEndfCont(); 
	      Process(secContNL);
	      sec->Add(secContNL);
	      for(Int_t l=0; l<secContNL->GetN1(); ++l) {
		TNudyEndfList *secList = new TNudyEndfList();
		Process(secList);
		sec->Add(secList); 
	      }
	    }
	    break;
	  case 7: // R-Matrix Limited Format
	    {
	      TNudyEndfCont *secCont=new TNudyEndfCont(); 
	      Process(secCont);
	      sec->Add(secCont);
	      
	      TNudyEndfCont *secContNJSX=new TNudyEndfCont(); 
	      Process(secContNJSX);
	      sec->Add(secContNJSX);
	      for(Int_t l=0; l<secContNJSX->GetL1(); ++l) {
		TNudyEndfList *secList = new TNudyEndfList();
		Process(secList);
		sec->Add(secList);
	      }
	      TNudyEndfList *secList = new TNudyEndfList();
	      Process(secList);
	      sec->Add(secList);
	    }
	    break;
	  default:
	    Fatal("ProcessF32(TNudyEndfSec*)","Unknown values of LRF %d",iLRF);
	    break;
	  }
	}
	  break;
	case 2: { // LCOMP == 2
	  // Resolved resonance compact covariance format
	  switch(iLRF) {
	  case 1: // SLBW
	  case 2: // MLBW
	  case 3: // Reich-Moore
	    {
	      TNudyEndfList *secListNRSA = new TNudyEndfList();
	      Process(secListNRSA);
	      sec->Add(secListNRSA);
	      TNudyEndfCont *secContNNN = new TNudyEndfCont();
	      Process(secContNNN);
	      sec->Add(secContNNN);
	      Int_t iNM = secContNNN->GetN1();
	      Int_t iNDIGIT = secContNNN->GetL1();
	      for(Int_t i1 = 0; i1 < iNM; i1++){
		TNudyEndfINTG *secINTG = new TNudyEndfINTG(18,iNDIGIT);
		Process(secINTG);
		sec->Add(secINTG);
	      }
	    }
	    break;
	  case 7:
	    {
	      TNudyEndfList *secListNPP = new TNudyEndfList();
	      Process(secListNPP);
	      sec->Add(secListNPP);
	      TNudyEndfList *secListNCH = new TNudyEndfList();
	      Process(secListNCH);
	      sec->Add(secListNCH);
	      TNudyEndfList *secListNX = new TNudyEndfList();
	      Process(secListNX);
	      sec->Add(secListNX);
	      
	      TNudyEndfCont *secContNNNB = new TNudyEndfCont();
	      Process(secContNNNB);
	      sec->Add(secContNNNB);
	      Int_t iNM = secContNNNB->GetN1();
	      Int_t iNDIGIT = secContNNNB->GetL1();

	      for(Int_t i2 = 0; i2 < iNM; i2++){

		TNudyEndfINTG *secINTG = new TNudyEndfINTG(18,iNDIGIT);
		Process(secINTG);
		sec->Add(secINTG);
	      }
	    }
	    break;
	  default:
	    Fatal("ProcessF32(TNudyEndfSec*)","Unknown values of LRF %d",iLRF);
	    break;
	  }
	}
	  break;
	default:
	  Fatal("ProcessF32(TNudyEndfSec*)","Unknown values of LCOMP %d",iLCOMP);
	  break;
	}
      } 
	break;
      case 2: { // LRU == 2
	TNudyEndfCont *secCont=new TNudyEndfCont(); 
	Process(secCont);
	sec->Add(secCont);
	secCont->Print(" ");///////////////////////////////////////////////test
	for(Int_t k=0; k<secCont->GetN1(); ++k) {
	  TNudyEndfList *secList = new TNudyEndfList();
	  Process(secList);
	  sec->Add(secList);
	}
	TNudyEndfList *secList = new TNudyEndfList();
	Process(secList);
	sec->Add(secList);
      }
	break;
      default:
	Fatal("ProcessF32(TNudyEndfSec*)","Unknown values of LRU %d",iLRU);
	break;
      }
    }
  }
  GetSEND(mtf);
}

//_______________________________________________________________________________
void TNudyENDF::ProcessF33(TNudyEndfSec *sec)
{
  Double_t c[2];
  Int_t  nl[4];
  Int_t mtf[4];
  Int_t &curMF  = mtf[1];
  Int_t &iNL = nl[3];

  GetCONT(c,nl,mtf);
  if(curMF!=33) 
    Fatal("ProcessF33(TNudyEndfSec*)","File %d should not be processed here",curMF);

  for(Int_t i=0; i<iNL; ++i) {
    TNudyEndfCont *secCont=new TNudyEndfCont(); 
    Process(secCont);
    sec->Add(secCont);
    for(Int_t j=0; j<secCont->GetN1(); ++j) {
      // NC-Type records
      TNudyEndfCont *secContNC=new TNudyEndfCont(); 
      Process(secContNC);
      sec->Add(secContNC);

      TNudyEndfList *secList = new TNudyEndfList();
      Process(secList);
      sec->Add(secList);
    }
    for(Int_t k=0; k<secCont->GetN2(); ++k) {
      // NI-Type records
      TNudyEndfList *secList = new TNudyEndfList();
      Process(secList);
      sec->Add(secList);
    }
  }
  GetSEND(mtf);
}

//_______________________________________________________________________________
void TNudyENDF::ProcessF34(TNudyEndfSec *sec)
{
  Double_t c[2];
  Int_t  nl[4];
  Int_t mtf[4];
  Int_t &curMF  = mtf[1];
  Int_t iNL = -1;
  Int_t iNL1 = -1;

  GetCONT(c,nl,mtf);
  if(curMF!=34) 
    Fatal("ProcessF34(TNudyEndfSec*)","File %d should not be processed here",curMF);

  for(Int_t i=0; i<sec->GetN2(); ++i) {
    TNudyEndfCont *secContMT=new TNudyEndfCont(); 
    Process(secContMT);
    sec->Add(secContMT);
    iNL = secContMT->GetN1();
    iNL1 = secContMT->GetN2();
    for(Int_t j=0; j<iNL; ++j) {
      for(Int_t k=j; k<iNL1; ++k){
        TNudyEndfCont *secCont=new TNudyEndfCont(); 
        Process(secCont);
        sec->Add(secCont);
        for(Int_t l=0; l<secCont->GetN2(); ++l) {
          TNudyEndfList *secList = new TNudyEndfList();
          Process(secList);
          sec->Add(secList);
        }
      }
    }
  }
  GetSEND(mtf);
}

//_______________________________________________________________________________
void TNudyENDF::ProcessF35(TNudyEndfSec *sec)
{
  Double_t c[2];
  Int_t  nl[4];
  Int_t mtf[4];
  Int_t &curMF  = mtf[1];

  GetCONT(c,nl,mtf);
  if(curMF!=35) 
    Fatal("ProcessF35(TNudyEndfSec*)","File %d should not be processed here",curMF);

  for(Int_t i=0; i<sec->GetN1(); ++i){
    TNudyEndfList *secList = new TNudyEndfList();
    Process(secList);
    sec->Add(secList);
  }  
  GetSEND(mtf);
}

//_______________________________________________________________________________
void TNudyENDF::ProcessF40(TNudyEndfSec *sec)
{
  Double_t c[2];
  Int_t  nl[4];
  Int_t mtf[4];
  Int_t &curMF  = mtf[1];

  GetCONT(c,nl,mtf);
  if(curMF!=40) 
    Fatal("ProcessF40(TNudyEndfSec*)","File %d should not be processed here",curMF);

  for(Int_t i=0; i<sec->GetN1(); ++i){
    TNudyEndfCont *secContNL=new TNudyEndfCont(); 
    Process(secContNL);
    sec->Add(secContNL);
    for(Int_t j=0; j<secContNL->GetN2(); ++j){
      TNudyEndfCont *secContSub=new TNudyEndfCont(); 
      Process(secContSub);
      sec->Add(secContSub);
      for(Int_t k=0; k<secContSub->GetN1(); ++k){
        // NC-Type records
        TNudyEndfCont *secCont=new TNudyEndfCont(); 
        Process(secCont);
        sec->Add(secCont);

        TNudyEndfList *secList = new TNudyEndfList();
        Process(secList);
        sec->Add(secList);
      }
      for(Int_t l=0; l<secContSub->GetN2(); ++l){
        // NI-Type records
        TNudyEndfList *secList = new TNudyEndfList();
        Process(secList);
        sec->Add(secList);
      }
    }
  }  
  GetSEND(mtf);
}

//_______________________________________________________________________________
void TNudyENDF::ToEndSec()
{
  Int_t mtf[4];
  Int_t &curMT  = mtf[2];

  while(!fENDF.eof()) {
    fENDF.getline(fLine,LINLEN);
    GetMTF(mtf);
    if(!curMT) 
      break;
  }
  // Get the Section END record and check it
  CheckSEND(mtf);
}

//_______________________________________________________________________________
void TNudyENDF::Process(TNudyEndfCont *secCont)
{
  Double_t c[2];
  Int_t nl[4];
  Int_t mtf[4];

  fENDF.getline(fLine,LINLEN);
  GetCONT(c,nl,mtf);
  secCont->SetCont(c[0],c[1],nl[0],nl[1],nl[2],nl[3]);

}


//_______________________________________________________________________________
void TNudyENDF::Process(TNudyEndfList *secList)
{
  Double_t c[2];
  Int_t nl[4];
  Int_t mtf[4];
  Double_t list[6];
  Int_t np=0;

  fENDF.getline(fLine,LINLEN);
  
  GetCONT(c,nl,mtf);

  secList->SetCont(c[0],c[1],nl[0],nl[1],nl[2],nl[3]);

  for(Int_t i=0; i<(secList->GetNPL()-1)/6+1; ++i) {
    fENDF.getline(fLine,LINLEN);
    GetFloat(list);

    for(Int_t j=0; j<6; ++j) {
      secList->SetLIST(list[j],np);
      if(++np >= secList->GetNPL()) break;
    }
  }
}

//_______________________________________________________________________________
void TNudyENDF::Process(TNudyEndfTab1 *secTab1)
{
  Double_t c[2];
  Int_t nl[4];
  Int_t mtf[4];
  Int_t nbt[6];
  Double_t xy[6];
  Int_t nr=0;
  Int_t np=0;
  fENDF.getline(fLine,LINLEN);
  GetCONT(c,nl,mtf);
  secTab1->SetCont(c[0],c[1],nl[0],nl[1],nl[2],nl[3]);
  if(mtf[0]==128 && mtf[1]==3 && mtf[2]==16) {
  }
  if(secTab1->GetNR() != 0)  {
    for(Int_t i=0; i<(secTab1->GetNR()-1)/3+1; ++i) {
      fENDF.getline(fLine,LINLEN);
      GetInt(nbt);
      for(Int_t j=0; j<6; j+=2) {
	secTab1->SetNBT(nbt[j],nr);
	secTab1->SetINT(nbt[j+1],nr);
	if(++nr >= secTab1->GetNR()) break;
      }
    }
  }
  if(secTab1->GetNP() != 0){
    for(Int_t i=0; i<(secTab1->GetNP()-1)/3+1; ++i) {
      fENDF.getline(fLine,LINLEN);
      GetFloat(xy);
      for(Int_t j=0; j<3; ++j) {
	secTab1->SetX(xy[2*j],np);
	secTab1->SetY(xy[2*j+1],np);
	if(++np >= secTab1->GetNP()) break;
      }
    }
  }
}

//_______________________________________________________________________________
void TNudyENDF::Process(TNudyEndfTab2 *secTab2)
{
  Double_t c[2];
  Int_t nl[4];
  Int_t mtf[4];
  Int_t nbt[6];
  Int_t nr=0;

  fENDF.getline(fLine,LINLEN);
  GetCONT(c,nl,mtf);
  secTab2->SetCont(c[0],c[1],nl[0],nl[1],nl[2],nl[3]);

  for(Int_t i=0; i<(secTab2->GetNR()-1)/3+1; ++i) {
    fENDF.getline(fLine,LINLEN);
    GetInt(nbt);
    for(Int_t j=0; j<6; j+=2) {
      secTab2->SetNBT(nbt[j],nr);
      secTab2->SetINT(nbt[j+1],nr);
      if(++nr >= secTab2->GetNR()) break;
    }
  }
}

void TNudyENDF::Process(TNudyEndfINTG *secINTG)
{
  Int_t mtf[4];
  Int_t IJ[2];
  Int_t KIJ[18]={0};
  Int_t ndigit = secINTG->GetNdigit();
  fENDF.getline(fLine,LINLEN);     
  if(ndigit < 2 || ndigit > 6) {    //  ****   ABHIJIT:   After complete C++ does the code need Ndigit?
    Error("TNudy::Process(TNudyEndfINTG *)","Invalid value for NDIGIT proceeding with NDIGIT = 2");
    secINTG->SetNdigit(2);
    ndigit = secINTG->GetNdigit();
  }
 
  GetINTG(ndigit,IJ,KIJ,mtf);
  secINTG->SetKIJ(KIJ);
  secINTG->SetIJ(IJ);
}

//_______________________________________________________________________________
void TNudyENDF::GetSEND(const Int_t pmtf[3])
{
  fENDF.getline(fLine,LINLEN);
  CheckSEND(pmtf);
}

//_______________________________________________________________________________
void TNudyENDF::GetFEND(const Int_t pmtf[3])
{
  fENDF.getline(fLine,LINLEN);
  CheckFEND(pmtf);
}

//_______________________________________________________________________________
void TNudyENDF::GetMEND(const Int_t pmtf[3])
{
  fENDF.getline(fLine,LINLEN);
  CheckMEND(pmtf);
}

//_______________________________________________________________________________
void TNudyENDF::GetTEND()
{
  fENDF.getline(fLine,LINLEN);
  CheckTEND();
}

//_______________________________________________________________________________
void TNudyENDF::CheckSEND(const Int_t pmtf[3]) const
{
  Double_t c[2];
  Int_t nl[4];
  Int_t mtf[4];

   GetCONT(c,nl,mtf);
 
  if( c[0]||c[1]||
     nl[0]||nl[1]||nl[2]||nl[3]||
    mtf[0]!=pmtf[0]||mtf[1]!=pmtf[1]||mtf[2] )
    {
      std::cout << "Expecting SEND and found for MAT " << pmtf[0] 
	   << " File " << pmtf[1] 
	   << " Section " << pmtf[2] << std::endl;
      Fatal("CheckSEND","SEND not found\n");
    }
}

//_______________________________________________________________________________
void TNudyENDF::CheckFEND(const Int_t pmtf[3]) const
{
  Double_t c[2];
  Int_t nl[4];
  Int_t mtf[4];

  GetCONT(c,nl,mtf);

  if(c[0]||c[1]||
       nl[0]||nl[1]||nl[2]||nl[3]||
     mtf[0]!=pmtf[0]||mtf[1]!=pmtf[1]||mtf[2]) 
    {
      std::cout << "Expecting FEND and found for MAT " << pmtf[0] 
	   << " File " << pmtf[1] << std::endl;
      Fatal("CheckFEND","FEND not found\n");
    }
}

//_______________________________________________________________________________
void TNudyENDF::CheckMEND(const Int_t pmtf[3]) const
{
  Double_t c[2];
  Int_t nl[4];
  Int_t mtf[4];

  GetCONT(c,nl,mtf);
  if(c[0]||c[1]||
     nl[0]||nl[1]||nl[2]||nl[3]||
     mtf[0]||mtf[1]||mtf[2]) 
    {
      std::cout << "Expecting MEND and found for MAT " << pmtf[0] << std::endl;
      Fatal("CheckMEND","MEND not found\n");
    }
}

//_______________________________________________________________________________
void TNudyENDF::CheckTEND() const
{
  Double_t c[2];
  Int_t nl[4];
  Int_t mtf[4];

  GetCONT(c,nl,mtf);
  if(c[0]||c[1]||
     nl[0]||nl[1]||nl[2]||nl[3]||
     (mtf[0]!=-1)||mtf[1]||mtf[2]) 
    {
      std::cout << "Expecting TEND and found" << std::endl;
      Fatal("CheckTEND","TEND not found\n");
    }
}

//_______________________________________________________________________________
void TNudyENDF::DumpENDF(Int_t flags=1)
{
  fTape->DumpENDF(flags);
}

//_______________________________________________________________________________
TNudyENDF::~TNudyENDF()
{
  //  printf("Endf destroyed");
  SafeDelete(fRENDF);
  SafeDelete(fTape);
}
