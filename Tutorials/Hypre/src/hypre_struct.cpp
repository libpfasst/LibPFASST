#include "hypre_struct.hpp"

void HypreStruct::SetComm(MPI_Comm in_comm)
{
   comm = in_comm;
   MPI_Comm_rank(comm, &myid);
   MPI_Comm_size(comm, &num_procs);
}

int HypreStruct::GetNumRows(void)
{
   return nrows;
}

void HypreStruct::SetDim(int in_dim)
{
   dim = in_dim;
}

void HypreStruct::InitGrid(int num_grid_points, int in_nrows, int *extents)
{
   ilower = (int *)malloc(2 * sizeof(int));
   iupper = (int *)malloc(2 * sizeof(int));
   stencil_indices = (int *)malloc(5 * sizeof(int));
   Lx = 1.0; Ly = 1.0;

   if (in_nrows < 0){
      n = num_grid_points;
      N  = sqrt(num_procs);
      h  = 1.0 / (N*n+1); /* note that when calculating h we must
                             remember to count the boundary nodes */
      h2 = h*h;
      pj = myid / N;
      pi = myid - pj*N;
    
      /* Figure out the extents of each processor's piece of the grid. */
      ilower[0] = pi*n;
      ilower[1] = pj*n;
    
      iupper[0] = ilower[0] + n-1;
      iupper[1] = ilower[1] + n-1;
    
      nrows = n*n;
    
      bc_nentries = 1;
      bc_nnz  = bc_nentries*n; /*  number of stencil entries times the length
                                     of one side of my grid box */
   }
   else {
      nrows = in_nrows;
      ilower[0] = extents[0];
      ilower[1] = extents[1];
      iupper[0] = extents[2];
      iupper[1] = extents[3];
      n = sqrt(nrows);
      N  = sqrt(num_procs);
      h  = 1.0 / (N*n+1);
   }

   nentries = 5;
   nnz = nentries*nrows;

   /* Set up a grid */
   {
      /* Create an empty grid object */
      HYPRE_StructGridCreate(comm, dim, &grid);

      /* Add a new box to the grid */
      HYPRE_StructGridSetExtents(grid, ilower, iupper);

      /* This is a collective call finalizing the grid assembly.
         The grid is now ``ready to be used'' */
      HYPRE_StructGridAssemble(grid);
   }

   /* Define the discretization stencil */
   {
      /* Create an empty stencil object */
      HYPRE_StructStencilCreate(2, nentries, &stencil);

      /* Define the geometry of the stencil */
      {
         for (int entry = 0; entry < nentries; entry++){
            HYPRE_StructStencilSetElement(stencil, entry, offsets_2D[entry]);
         }
      }
   }

 
   coords_x = (double *)malloc(nrows * sizeof(double));
   coords_y = (double *)malloc(nrows * sizeof(double));

   //int k = 0;
   //for (int i = ilower[0]; i <= iupper[0]; i++){
   //   for (int j = ilower[1]; j <= iupper[1]; j++){
   //      coords_x[k] = h * (i+1);
   //      coords_y[k] = h * (j+1);
   //      k++;
   //   }
   //}

   int k = 0;
   for(int j = 0; j < n; j++){
      for (int i = 0; i < n; i++){
         coords_x[k] = (ilower[0]+i) * h + h;
         coords_y[k] = (ilower[1]+j) * h + h;
         k++;
      }
   }
}

double HeatEquExactSin(double t, double coord, double nu, double kfreq, double L)
{
   double u;
   double omega = kfreq * M_PI / L;
   u = sin(omega * coord) * exp(-nu * pow(omega, 2.0) * t);
   return u;
}

double *HypreStruct::HeatEquTrueSol(double t, int P, int Q, double init_cond)
{
   double *u = (double *)calloc(nrows, sizeof(double));
  // for (int p = 1; p <= P; p++){
  //    for (int q = 1; q <= Q; q++){
  //       int k = 0;
  //       for (int i = ilower[0]; i <= iupper[0]; i++){
  //          for (int j = ilower[1]; j <= iupper[1]; j++){
  //             double a = q * M_PI / Lx, b = p * M_PI / Ly;
  //             double c = (2 / sqrt(Lx * Ly)) * init_cond * (cos(a*Lx) - 1)*(cos(b*Ly) - 1) / (a*b);
  //             u[k] = u[k] + c * sin(a * coords_x[k]) * sin(b * coords_y[k]) * exp(-pow(M_PI, 2.0) * (pow(q/Lx, 2.0) + pow(p/Ly, 2.0)) * t);
  //             k++;
  //          }
  //       }
  //    }
  // }
  // int k = 0;
  // for (int i = ilower[0]; i <= iupper[0]; i++){
  //    for (int j = ilower[1]; j <= iupper[1]; j++){
  //       u[k] = 2 * u[k] / sqrt(Lx * Ly);
  //       k++;
  //    }
  // }
   int k = 0;
   for (int i = ilower[0]; i <= iupper[0]; i++){
      for (int j = ilower[1]; j <= iupper[1]; j++){
         u[k] = HeatEquExactSin(t, coords_x[k], 1.0, 1.0, Lx) * HeatEquExactSin(t, coords_y[k], 1.0, 1.0, Ly);
         k++;
      }
   }
   return u;
}
