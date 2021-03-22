#ifndef CENCAP_H
#define CENCAP_H

#include "cscalar.h"

/*
 *  Prototypes for "Encap" functions called from Fortran to manipulate Scalar objects
 */

extern "C"
{
   void ScalarCreate(Scalar *x);
   void ScalarDestroy(Scalar *x);
   void ScalarSetval(Scalar *x, double y);
   void ScalarCopy(Scalar *x, Scalar *y);
   void ScalarPack(Scalar *x, double *y);
   void ScalarUnpack(Scalar *x, double y);
   void ScalarNorm(Scalar *x, double *norm);
   void ScalarAxpy(Scalar *y, double a, Scalar *x);
   void ScalarPrint(Scalar *x);
}

#endif

