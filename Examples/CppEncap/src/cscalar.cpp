#include "cscalar.hpp"
#include <iostream>
#include <cstdlib>

Scalar::Scalar(void)
{

}

Scalar::~Scalar(void)
{

}

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
