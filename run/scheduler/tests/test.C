void test()
{
  const Int_t mask = 0xFF;
  Int_t arr[256];
  for (UInt_t i = 0; i < 256; i++) {
    arr[i] = 8;
    for (UInt_t j = 8; j > 0; j--) {
      if (((1 << (j - 1)) & i) == 0) {
        arr[i] = j - 1;
        break;
      }
    }
  }
  for (UInt_t i = 0; i < 16; i++) {
    for (UInt_t j = 0; j < 16; j++) {
      UInt_t tt = 16 * i + j;
      printf("%d,", arr[tt]);
    }
    printf("\n");
  }
}
