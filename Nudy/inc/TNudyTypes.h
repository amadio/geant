#ifndef ROOT_TNudyTypes
#define ROOT_TNudyTypes

enum Reaction_t { kNoReaction = 0, kN_TOTAL = 1, kZ_Z0 = 2, kZ_NONELAS = 3, kZ_N = 4, kZ_ANY = 5, kZ_CONTIN = 10, kZ_2ND = 11, kZ_2N = 16, kZ_3N = 17, kZ_FISSION = 18,
		  kN_F = 19, kN_NF = 20, kN_2NF = 21, kZ_NALPHA = 22, kN_N3ALPHA = 23, kZ_2NALPHA = 24, kZ_3NALPHA = 25, kN_ABS = 27, kZ_NP = 28, kZ_N2ALPHA = 29,
		  kZ_2N2ALPHA = 30, kZ_ND = 32, kZ_NT = 33, kZ_N3HE = 34, kZ_ND2ALPHA = 35, kZ_NT2ALPHA = 36, kZ_4N = 37, kN_3NF = 38, kZ_2NP = 41, kZ_3NP = 42,
		  kZ_N2P = 44, kZ_NPALPHA = 45, kY_N0 = 50, kZ_N1 = 51, kZ_N2 = 52, kZ_N40 = 90, kZ_NC = 91, kN_DISAP = 101, kZ_GAMMA = 102, kZ_P = 103,
		  kZ_D = 104, kZ_T = 105, kZ_3HE = 106, kZ_ALPHA = 107, kZ_2ALPHA = 108, kZ_3ALPHA = 109, kZ_2P = 111, kZ_PALPHA = 112, kZ_T2ALPHA = 113,
		  kZ_D2ALPHA = 114, kZ_PD = 115, kZ_PT = 116, kZ_DALPHA = 117, kN_RES = 151, kZ_XN = 201, kZ_XGAMMA = 202, kZ_XP = 203, kZ_XD = 204, kZ_XT = 205,
		  kZ_X3HE = 206, kZ_XALPHA = 207, kZ_XPIPOSITIVE = 208, kZ_XPI0 = 209, kZ_XPINEGATIVE = 210, kZ_XUPOSITIVE = 211, kZ_XUNEGATIVE = 212,
		  kZ_XKPOSITIVE = 213, kZ_XK0LONG = 214, kZ_XK0SHORT = 215, kZ_XKNEGATIVE = 216, kZ_XPNEGATIVE = 217, kZ_XNNEGATIVE = 218, kK = 534, kL1 = 535,
		  kL2 = 536, kL3 = 537, kM1 = 538, kM2 = 539, kM3 = 540, kM4 = 541, kM5 = 542, kN1 = 543, kN2 = 544, kN3 = 545, kN4 = 546, kN5 = 547, kN6 = 548,
		  kN7 = 549, kO1 = 550, kO2 = 551, kO3 = 552, kO4 = 553, kO5 = 554, kO6 = 555, kO7 = 556, kO8 = 557, kO9 = 558, kP1 = 559, kP2 = 560, kP3 = 561,
		  kP4 = 562, kP5 = 563, kP6 = 564, kP7 = 565, kP8 = 566, kP9 = 567, kP10 = 568, kP11 = 569, kQ1 = 570, kQ2 = 571, kQ3 = 572, kZ_P0 = 600,
		  kZ_P1 = 601, kZ_P2 = 602, kZ_P3 = 603, kZ_P4 = 604, kZ_PC = 649, kZ_D0 = 650, kZ_D1 = 651, kZ_D2 = 652, kZ_DC = 699, kZ_T0 = 700, kZ_T1 = 701,
		  kZ_T2 = 702, kZ_TC = 749, kN_3HE0 = 750, kN_3HE1 = 751, kN_3HEC = 799, kZ_ALPHA0 = 800, kZ_ALPHA1 = 801, kZ_ALPHAC = 849, kZ_2N0 = 875,
		  kZ_2N1 = 876, kZ_2NC = 891 };

enum FileData_t {
	kGen_Info=1, kResonance_Param=2, kReac_XSect=3, kAng_Dist=4, kEnergy_Dist=5,
	kPE_AD=6, kTNSCLD=7, kDFPY=8, kMRP=9, kPXSectRN=10, kGCPhotonProd=11, kPPYD=12,
	kPPXSect=13, kPAng_Dist=14, kContPES=15, kPhotonInt=23, kSDPEAD=26, kAFFSF=27,
	kARD=28, kCMP=30, kCRF=31, kCRP=32, kNXSect=33, kCAng_Dist=34, kCE_Dist=35,
	kCRNP=40
};

enum AliasDist_t {
  kOriginal=-1, kBuilt=1
};

enum IntScheme_t {
  kLinear = 1, kDiscrete = 2, kUniform = 3
};

#define ERROR_MARGIN 1e-5

#define TNUDYALIAS_MULTITHREAD

#endif
