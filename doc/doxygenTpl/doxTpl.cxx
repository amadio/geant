#include "doxTpl.h"

/**
   @file doxTpl.cxx Implementation file for doxTpl

   The \@file tag should be the only Doxygen tag in the implementation file 
   since we would like to keep all the documentation in the same place for
   ease of access and formatting.
*/

//________________________________________________________________________________
doxTpl::doxTpl() {
   // Non doxygen comment that will not be stripped from the source
}

//________________________________________________________________________________
int doxTpl::oneMethod(int a) {
   // Another local comment
   int b=2*a;
   return b;
}

