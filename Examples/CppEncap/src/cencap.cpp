#include "cencap.hpp"
#include <new>

/*
 *  "Encap" functions called from Fortran to manipulate Scalar objects
 */

extern "C"
{
   void ScalarCreate(Scalar **x)
   {
      *x = new Scalar();
   }

   void ScalarDestroy(Scalar *x)
   {
      if (x != nullptr){
         delete x;
      }
   }

   void ScalarSetVal(Scalar *x, double y)
   {
      x->DataSetVal(y);
   }

   void ScalarCopy(Scalar *dest, Scalar *src)
   {
      dest->DataSetVal(src->DataGet());
   }

   double ScalarPack(Scalar *x)
   {
      return x->DataGet();
   }

   void ScalarUnpack(Scalar *x, double y)
   {
      x->DataSetVal(y);
   }

   double ScalarNorm(Scalar *x)
   {
      return x->DataNorm();
   }

   void ScalarAxpy(Scalar *y, double a, Scalar *x)
   {
      y->DataAxpy(a, x->DataGet());
   }

   void ScalarPrint(Scalar *x)
   {
      x->DataPrint();
   }

   double ScalarGetVal(Scalar *x)
   {
      return x->DataGet();
   }
}

