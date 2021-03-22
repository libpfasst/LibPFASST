#ifndef CENCAP_HPP
#define CENCAP_HPP

#include "cscalar.hpp"

/*
 *  Prototypes for "Encap" functions called from Fortran to manipulate Scalar objects
 */

extern "C"
{
   void ScalarCreate(Scalar **x);
   void ScalarDestroy(Scalar *x);
   void ScalarSetVal(Scalar *x, double y);
   void ScalarCopy(Scalar *dest, Scalar *src);
   double ScalarPack(Scalar *x);
   void ScalarUnpack(Scalar *x, double y);
   double ScalarNorm(Scalar *x);
   void ScalarAxpy(Scalar *y, double a, Scalar *x);
   void ScalarPrint(Scalar *x);
   double ScalarGetVal(Scalar *x);
}

#endif
