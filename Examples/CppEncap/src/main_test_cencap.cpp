#include "cscalar.hpp"
#include "cencap.hpp"

int main(void)
{
   Scalar *s, *x;
   ScalarCreate(&s);
   ScalarCreate(&x);

   double s_val = -1.6;
   double a = 2.8, x_val = -.43;

   ScalarSetVal(s, s_val);
   cout << "True value is " << s_val << ".  ScalarGetVal() returns " << ScalarGetVal(s) << ".\n";
   double norm = ScalarNorm(s);
   cout << "True norm is " << abs(s_val) << ". ScalarNorm() returns " << norm << ".\n";
   double axpy = s_val + a*x_val;
   ScalarSetVal(x, x_val);
   ScalarAxpy(s, a, x);
   cout << "True axpy is " << axpy << ".  ScalarGetVal() returns " << ScalarGetVal(s) << " after ScalarAxpy() called.\n"; 

   ScalarDestroy(s);
   ScalarDestroy(x);

   return 0;
}
