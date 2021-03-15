#include "hypre_vector.hpp"
#include "hypre_encap.hpp"

int main(void)
{
   HypreVector *x, *y;
   int num_grid_points = 3;
   
   MPI_Init(NULL, NULL);

   double a = 2.8;
   x = new HypreVector(num_grid_points, 0.0, MPI_COMM_WORLD);
   y = new HypreVector(num_grid_points, 0.0, MPI_COMM_WORLD);
   
   x->SetVal(1.0);
   y->SetVal(2.0);
   y->Axpy(a, x);

   delete x;
   delete y;

   HypreVectorCreate(&x, num_grid_points, MPI_COMM_WORLD);
   HypreVectorCreate(&y, num_grid_points, MPI_COMM_WORLD);

   HypreVectorSetVal(x, 1.0);
   HypreVectorSetVal(y, 2.0);
   
   HypreVectorAxpy(y, a, x);

   printf("\n%f %f\n", HypreVectorNorm(x), HypreVectorNorm(y));

   HypreVectorPrint(y);
   HypreVectorCopy(x, y);
   HypreVectorPrint(y);
   double *values = (double *)malloc(y->GetNumRows() * sizeof(double));
   for (int i = 0; i < y->GetNumRows(); i++){
      values[i] = a;
   }
   HypreVectorUnpack(y, values);
   HypreVectorPrint(y);

   HypreVectorSetVal(y, 3.0);
   double *z = HypreVectorPack(y);
   printf("\n");
   for (int i = 0; i < y->GetNumRows(); i++){
      printf("%f\n", z[i]);
   }

   HypreVectorDestroy(x);
   HypreVectorDestroy(y);

   MPI_Finalize();
   return 0;
}
