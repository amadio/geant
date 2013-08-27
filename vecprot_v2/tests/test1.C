#include "TBits.h"
#include "TRandom.h"

void test1()
{
   const UInt_t size = 100000;
   TBits b(size);
   Bool_t success = kTRUE;
   UInt_t i, j;
   printf("Testing TBits::LastSetBit() ...\n");
   for (i=0; i<size; i++) {
      b.SetBitNumber(i);
      if ((b.LastSetBit() != i) || (b.LastSetBit(i+1) != i)) {
         printf("Error for bit %d: lastset=%d lastset(i+1)=%d\n", i, b.LastSetBit(), b.LastSetBit(i+1));
         success = kFALSE;
      }   
      b.SetBitNumber(i, kFALSE);
   }
   if (success) printf("Test 1: Set/Find SUCCESS\n");
   else         printf("Test 1: Set/Find FAILED\n");
   
   for (i=0; i<size; i++) {
      j = size*gRandom->Rndm();
      b.SetBitNumber(j);
   }
   Long64_t checksum1 = 0;
   j = b.FirstSetBit();
   while (j<size) {
      checksum1 += j;
      j = b.FirstSetBit(j+1);
   }
   Long64_t checksum2 = 0;
   j = b.LastSetBit();
   while (j<size) {
      checksum2 += j;
      if (j==0) break;
      j = b.LastSetBit(j-1);
   }
   printf("Test 2: checksum1=%lld  checksum2=%lld ... ", checksum1, checksum2);
   if (checksum1==checksum2) printf("SUCCESS\n");
   else         printf("FAILED\n");   

   printf("Testing TBits::LastNullBit() ...\n");
   for (i=0; i<size; i++) b.SetBitNumber(i);
   for (i=0; i<size; i++) {
      b.SetBitNumber(i, false);
      if ((b.LastNullBit() != i) || (b.LastNullBit(i+1) != i)) {
         printf("Error for bit %d: lastnull=%d lastnull(i+1)=%d\n", i, b.LastNullBit(),b.LastNullBit(i+1));
         success = kFALSE;
      }   
      b.SetBitNumber(i);
   }
   if (success) printf("Test 2: Reset/Find SUCCESS\n");
   else         printf("Test 2: Reset/Find FAILED\n");
   
   for (i=0; i<size; i++) {
      j = size*gRandom->Rndm();
      b.SetBitNumber(j, false);
   }
   Long64_t checksum3 = 0;
   j = b.FirstNullBit();
   while (j<size) {
      checksum3 += j;
      j = b.FirstNullBit(j+1);
   }
   Long64_t checksum4 = 0;
   j = b.LastNullBit();
   while (j<size) {
      checksum4 += j;
      if (j==0) break;
      j = b.LastNullBit(j-1);
   }
   printf("Test 2: checksum3=%lld  checksum4=%lld ... ", checksum3, checksum4);
   if (checksum3==checksum4) printf("SUCCESS\n");
   else         printf("FAILED\n");   
}
