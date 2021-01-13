#include "hypre_vector.hpp"

HypreVector::HypreVector(int num_grid_points, double value, MPI_Comm in_comm, int in_dim, int in_nrows, int *extents)
{
   SetComm(in_comm);
   SetDim(in_dim);
   InitGrid(num_grid_points, in_nrows, extents);
   //printf("%d %d %d %d %d\n", nrows, ilower[0], ilower[1], iupper[0], iupper[1]);

   double *values = (double *)calloc(nrows, sizeof(double));

   HYPRE_StructVectorCreate(comm, grid, &v);
   HYPRE_StructVectorInitialize(v);

   for (int i = 0; i < nrows; i++){
      values[i] = value;
   }
   HYPRE_StructVectorSetBoxValues(v, ilower, iupper, values);

   HYPRE_StructVectorAssemble(v);
}

HypreVector::~HypreVector(void)
{
   HYPRE_StructVectorDestroy(v);
}

HYPRE_StructVector HypreVector::GetVector(void)
{
   return v;
}

double *HypreVector::GetBoxValues(void)
{
   double *values = (double *)calloc(nrows, sizeof(double));
   HYPRE_StructVectorGetBoxValues(v, ilower, iupper, values);
   return values;
}

void HypreVector::SetBoxValues(double *values)
{
   HYPRE_StructVectorSetBoxValues(v, ilower, iupper, values);
}

void HypreVector::SetVal(double val)
{
   HYPRE_StructVectorSetConstantValues(v, val);
}

void HypreVector::SetSinInitCond(void)
{
   double *values = (double *)calloc(nrows, sizeof(double));
   double kfreq = 1.0;
   double omega_x = kfreq * M_PI / Lx;
   double omega_y = kfreq * M_PI / Ly;



   if (dim == 1){
      int k = 0;
      for (int i = ilower[0]; i <= iupper[0]; i++){
         values[k] = sin(omega_x * coords_x[k]);
         k++;
      }
   }
   else {
      int k = 0;
      for (int i = ilower[0]; i <= iupper[0]; i++){
         for (int j = ilower[1]; j <= iupper[1]; j++){
            values[k] = sin(omega_x * coords_x[k]) * sin(omega_y * coords_y[k]);
            k++;
         }
      }
   }
   HYPRE_StructVectorSetBoxValues(v, ilower, iupper, values);
}

void HypreVector::Axpy(double alpha, HypreVector *x)
{
   hypre_StructAxpy(alpha, x->GetVector(), v);
}

double HypreVector::Norm(void)
{
   return sqrt(hypre_StructInnerProd(v, v));
}

void HypreVector::Print(void)
{
   double *values = GetBoxValues();
   if (myid == 0) printf("\n");
   for (int p = 0; p < num_procs; p++){
      if (p == myid){  
         for (int i = 0; i < nrows; i++){
            printf("%.16f\n", values[i]);
         }
      }
      MPI_Barrier(comm);
   }
   if (myid == 0) printf("\n");
  // HYPRE_StructVectorPrint("hypre_vector_data.txt", v, 1);
}
