#ifndef CSCALAR_HPP
#define CSCALAR_HPP

#include <iostream>
#include <cstdlib>

using namespace std;

class Scalar {
   protected:
      double data;
   public:
      Scalar(void);
      ~Scalar(void);

      double DataGet(void);
      void DataSetVal(double y);
      double DataNorm(void);
      void DataAxpy(double a, double x);
      void DataPrint(void);
};

#endif
