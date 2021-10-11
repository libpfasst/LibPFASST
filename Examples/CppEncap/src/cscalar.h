#ifndef CSCALAR_H
#define CSCALAR_H

#include <iostream>
#include <cstdlib>

using namespace std;

class Scalar {
   protected:
      double data;
   public:
      Scalar(){
         ;
      }

      ~Scalar(){
         ;
      }

      double DataGet(void);
      void DataSetVal(double y);
      double DataNorm(void);
      void DataAxpy(double a, double x);
      void DataPrint(void);
};

double Scalar::DataGet(void)
{
   return data;
}

void Scalar::DataSetVal(double y)
{
   data = y;
}

double Scalar::DataNorm(void)
{
   return abs(data);
}

void Scalar::DataAxpy(double a, double x)
{
   data += a * x;
}

void Scalar::DataPrint()
{
   cout << data << '\n';
}

#endif
