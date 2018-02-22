#include "TBits.h"
#include "TRandom.h"

void Compact(TBits& b, Int_t nactive);
void PrintBits(TBits& b, Int_t nactive);
Int_t Reshuffle(TBits& fSelected, Int_t nactive);

void testcompact()
{
// Test compacting a vector of bits
   const Int_t size = 100;
   TBits b(size);
   // Generate some holes
   Int_t nactive = (size+1)*gRandom->Rndm();
   Int_t maxholes = (nactive+1)*gRandom->Rndm();
   for (Int_t i=0; i<maxholes; i++) b.SetBitNumber(nactive*gRandom->Rndm());
   for (Int_t i=nactive; i<size; i++) b.SetBitNumber(i);
//   Compact(b, nactive);
   Int_t nselected = Reshuffle(b, nactive);
   printf("nselected=%d\n", nselected);
}

void Compact(TBits& fHoles, Int_t nactive)
{
// Same algorithm as in Track
   PrintBits(fHoles, nactive);
   if (nactive == 0) return;
   Int_t firsthole = fHoles.FirstSetBit();
   while (firsthole<nactive) {
      Int_t lastactive = fHoles.LastNullBit(nactive-1);
      if (lastactive < nactive) {
         nactive = lastactive+1;
         if (firsthole==nactive) {
            PrintBits(fHoles, nactive);
            return;
         }
      } else {
         // No active tracks left
         nactive = 0;
         PrintBits(fHoles, nactive);
         return;
      }
      // exchange tracks pointed by firsthole and lastactive
//      SwapTracks(firsthole, lastactive);
      fHoles.SetBitNumber(firsthole, false);
      fHoles.SetBitNumber(lastactive, true);
      firsthole = fHoles.FirstSetBit(firsthole+1);
      nactive--;
      PrintBits(fHoles, nactive);
   }
}

Int_t Reshuffle(TBits& fSelected, Int_t nactive)
{
   PrintBits(fSelected, nactive);
   if (nactive == 0) return 0;
   Int_t nselected = nactive;
   Int_t firsthole = fSelected.FirstNullBit();
   while (firsthole<nselected) {
      Int_t lastsel = fSelected.LastSetBit(nselected-1);
      if (lastsel >= nselected) {PrintBits(fSelected, nselected); return 0;}
      nselected = lastsel+1;
      if (firsthole==nselected) {PrintBits(fSelected, nselected);return nselected;}
      // exchange tracks pointed by firsthole and lastactive
//      SwapTracks(firsthole, lastsel);
      fSelected.SetBitNumber(firsthole, true);
      fSelected.SetBitNumber(lastsel, false);
      firsthole = fSelected.FirstNullBit(firsthole+1);
      nselected--;
      PrintBits(fSelected, nselected);
   }
   return nselected;
}

void PrintBits(TBits& b, Int_t nactive)
{
   Int_t nbits = b.GetNbits();
   for (Int_t i=0; i<nbits; i++) {
      Int_t bit = (Int_t)b.TestBitNumber(i);
      printf("%d ", bit);
      if (i==nactive-1) printf("|| ");
   }
   printf("\n");
}      
