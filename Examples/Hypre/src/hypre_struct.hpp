#ifndef HYPRE_STRUCT_HPP
#define HYPRE_STRUCT_HPP

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include "_hypre_struct_mv.h"
#include "_hypre_struct_ls.h"
#include "HYPRE_struct_mv.h"
#include "HYPRE_struct_ls.h"
#include <mpi.h>

using namespace std;

class HypreStruct {
   protected:
      ;

   public:
      /* TODO: make protected */
      HYPRE_StructGrid grid;
      HYPRE_StructStencil stencil;
      int offsets_1D[5][2] = {{0,0}, {-1,0}, {1,0}, {0,0}, {0,0}};
      int offsets_2D[5][2] = {{0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}};

      double Lx, Ly;
      double h, h2;
      int n, N;
      int pi, pj;
      int **p_map;
      int nentries, bc_nentries;
      int nnz, nrows;
      int bc_nnz, bc_nrows;
      int *ilower, *iupper, bc_ilower[2], bc_iupper[2];
      int *stencil_indices, bc_stencil_indices[1];
      double *coords_x, *coords_y;
      double *exact_sol;
      MPI_Comm comm;
      int myid, num_procs;
      int dim;
      /* end TODO */


      HypreStruct() {
         ;
      }
      ~HypreStruct() {
         ;
      }

      void SetComm(MPI_Comm target_comm);
      void SetDim(int in_dim);
      int GetNumRows(void);
      void SetInitCond(double val);
      void InitGrid(int num_grid_points, int in_nrows = -1, int *extents = NULL);
      double *HeatEquTrueSol(double t, int P, int Q, double init_cond);
};

#endif
