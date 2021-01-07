#include "hypre_solver.hpp"

HypreSolver::HypreSolver(MPI_Comm in_comm, int in_dim, int in_max_iter, int in_max_levels)
{
   //if (setup_done == 1) return;
   SetComm(in_comm);
   SetDim(in_dim);

   n_pre  = 1;
   n_post = 1;
   tol = 0.0;
   max_iter = in_max_iter;
   jacobi_weight = 1.0;
   relax_type = RELAX_RBGS;
   solver_type = SOLVER_PFMG;
   max_levels = in_max_levels;
   setup_done = 1;
}

HypreSolver::~HypreSolver(void)
{
   Cleanup();
}

HYPRE_StructSolver HypreSolver::GetStructSolver(void)
{
   return solver;
}

int HypreSolver::GetNumRowsLevel(int level_index)
{
   return nrows_l[level_index];
}

int HypreSolver::GetNumLevels(void)
{
   return num_levels;
}

void HypreSolver::SetupMatrix(int num_grid_points, int level_index, int spacial_coarsen_flag)
{
   if ((level_index > 0) && (A == NULL) && (spacial_coarsen_flag == 1)) return;

   double *values;

   InitGrid(num_grid_points);

   /* Set up a Struct Matrix */
   {
      int nvalues = nnz;

      /* Create an empty matrix object */
      HYPRE_StructMatrixCreate(comm, grid, stencil, &A);

      /* Indicate that the matrix coefficients are ready to be set */
      HYPRE_StructMatrixInitialize(A);

      values = (double *)calloc(nvalues, sizeof(double));

      for (int j = 0; j < nentries; j++){
         stencil_indices[j] = j;
      }

      /* Set the standard stencil at each grid point,
         we will fix the boundaries later */
      if (dim == 1){
          for (int i = 0; i < nvalues; i += nentries){
            values[i] = -2.0 / h;
            for (int j = 1; j < nentries; j++){
               values[i+j] = 1.0 / h;
            }
         }
      }
      else {
         for (int i = 0; i < nvalues; i += nentries){
            values[i] = -4.0 / h2;
            for (int j = 1; j < nentries; j++){
               values[i+j] = 1.0 / h2;
            }
         }
      }

      HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, nentries,
                                     stencil_indices, values);

      free(values);
   }

   /* Incorporate the zero boundary conditions: go along each edge of
      the domain and set the stencil entry that reaches to the boundary to
      zero.*/
   {
      int bc_nvalues = bc_nnz;

      values = (double*) calloc(bc_nvalues, sizeof(double));
      for (int j = 0; j < bc_nvalues; j++){
         values[j] = 0.0;
      }
      if (dim == 1){
         if (pi == 0){
            /* Left grid point */
            bc_ilower[0] = 0;
            bc_iupper[0] = 0;
            bc_ilower[1] = 0;
            bc_iupper[1] = 0;

            bc_stencil_indices[0] = 1;

            HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, bc_nentries,
                                           bc_stencil_indices, values);
         }

         if (pi == N-1){
            /* Right grid point */
            bc_ilower[0] = pi*n + n-1;
            bc_iupper[0] = bc_ilower[0];
            bc_ilower[1] = 0;
            bc_iupper[1] = 0;

            bc_stencil_indices[0] = 2;

            HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, bc_nentries,
                                           bc_stencil_indices, values);
         }
      }
      else {
         /* Recall: pi and pj describe position in the processor grid */
         if (pj == 0){
            /* Bottom row of grid points */
            bc_ilower[0] = pi*n;
            bc_ilower[1] = pj*n;

            bc_iupper[0] = bc_ilower[0] + n-1;
            bc_iupper[1] = bc_ilower[1];

            bc_stencil_indices[0] = 3;

            HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, bc_nentries,
                                           bc_stencil_indices, values);
         }

         if (pj == N-1){
            /* upper row of grid points */
            bc_ilower[0] = pi*n;
            bc_ilower[1] = pj*n + n-1;

            bc_iupper[0] = bc_ilower[0] + n-1;
            bc_iupper[1] = bc_ilower[1];

            bc_stencil_indices[0] = 4;

            HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, bc_nentries,
                                           bc_stencil_indices, values);
         }

         if (pi == 0){
            /* Left row of grid points */
            bc_ilower[0] = pi*n;
            bc_ilower[1] = pj*n;

            bc_iupper[0] = bc_ilower[0];
            bc_iupper[1] = bc_ilower[1] + n-1;

            bc_stencil_indices[0] = 1;

            HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, bc_nentries,
                                           bc_stencil_indices, values);
         }

         if (pi == N-1){
            /* Right row of grid points */
            bc_ilower[0] = pi*n + n-1;
            bc_ilower[1] = pj*n;

            bc_iupper[0] = bc_ilower[0];
            bc_iupper[1] = bc_ilower[1] + n-1;

            bc_stencil_indices[0] = 2;

            HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, bc_nentries,
                                           bc_stencil_indices, values);
         }
      }
      free(values);
   }

   HYPRE_StructMatrixAssemble(A);

   /* Set up Struct Vectors for b and x */
   {
      int nvalues = nrows;
      double *values = (double*) calloc(nvalues, sizeof(double));

      /* Create an empty vector object */
      HYPRE_StructVectorCreate(comm, grid, &b);
      HYPRE_StructVectorCreate(comm, grid, &x);

      /* Indicate that the vector coefficients are ready to be set */
      HYPRE_StructVectorInitialize(b);
      HYPRE_StructVectorInitialize(x);

     /* Set the values */
      for (int i = 0; i < nvalues; i ++){
         values[i] = 0.0;
      }
      HYPRE_StructVectorSetBoxValues(b, ilower, iupper, values);
      HYPRE_StructVectorSetBoxValues(x, ilower, iupper, values);

      free(values);
   }

   HYPRE_StructVectorAssemble(b);
   HYPRE_StructVectorAssemble(x);

   if (level_index == 0){
      if (spacial_coarsen_flag == 1){
         SetupStructPFMGSolver(&pfmg_level_data); 
         HYPRE_StructPFMGSetMaxLevels(pfmg_level_data, max_levels);
         HYPRE_StructPFMGSetup(pfmg_level_data, A, b, x);
         pfmg_data = (hypre_PFMGData *)pfmg_level_data;
         A_l = pfmg_data->A_l;   
         RT_l = pfmg_data->RT_l;
         P_l = pfmg_data->P_l;
         x_l = pfmg_data->x_l;
         b_l = pfmg_data->b_l;
         restrict_data_l = pfmg_data->restrict_data_l;
         interp_data_l = pfmg_data->interp_data_l;
         num_levels = pfmg_data->num_levels;
      }
      else {
         num_levels = max_levels;      
      }
      nrows_l = (int *)malloc(num_levels * sizeof(int));
      ilower_l = (int **)malloc(num_levels * sizeof(int *));
      iupper_l = (int **)malloc(num_levels * sizeof(int *));
      for (int level = 0; level < num_levels; level++){
         ilower_l[level] = (int *)malloc(2 * sizeof(int));
         iupper_l[level] = (int *)malloc(2 * sizeof(int)); 
         HYPRE_StructMatrix AA;
         if (spacial_coarsen_flag == 1){
            AA = A_l[level];
         }
         else {
            AA = A;
         }
         nrows_l[level] = AA->grid->local_size;
         if (nrows_l[level] > 0){
            ilower_l[level] = AA->grid->boxes->boxes[0].imin; 
            iupper_l[level] = AA->grid->boxes->boxes[0].imax;
         }
         else {
            ilower_l[level][0] = 0; 
            ilower_l[level][1] = 0;
            iupper_l[level][0] = -1;
            iupper_l[level][1] = -1;
         }
      }
   }
}

void HypreSolver::FEval(double *y, double t, int level_index, double **f)
{
   double *values = (double *)calloc(nrows, sizeof(double));

   int num_iterations;
   double final_res_norm;

   for (int i = 0; i < nrows; i++){
      values[i] = 0.0;
   }
   HYPRE_StructVectorSetBoxValues(b, ilower, iupper, values);

   for (int i = 0; i < nrows; i++){
      values[i] = y[i];
   }
   HYPRE_StructVectorSetBoxValues(x, ilower, iupper, values);

   HYPRE_StructVectorAssemble(b);
   HYPRE_StructVectorAssemble(x);

   double alpha = 1.0, beta = 0.0;
   hypre_StructMatvec(alpha, A, x, beta, b);

   HYPRE_StructVectorGetBoxValues(b, ilower, iupper, values);
   for (int i = 0; i < nrows; i++){
      (*f)[i] = values[i];
   }

   free(values);
}

void HypreSolver::FComp(double **y, double t, double dtq, double *rhs, int level_index, double **f)
{
   int num_iterations;
   double final_res_norm;

   double *values = (double *)calloc(nnz, sizeof(double));
   double *space_values = (double *)calloc(nnz, sizeof(double));

   for (int i = 0; i < nrows; i++){
      values[i] = rhs[i];
   }
   HYPRE_StructVectorSetBoxValues(b, ilower, iupper, values);

   for (int i = 0; i < nrows; i++){
      values[i] = 0.0;
   }
   HYPRE_StructVectorSetBoxValues(x, ilower, iupper, values);

   HYPRE_StructVectorAssemble(b);
   HYPRE_StructVectorAssemble(x);

   HYPRE_StructMatrixGetBoxValues(A, ilower, iupper, nentries,
                                  stencil_indices, space_values);

   /* Update the matrix entries using temporal information */
  // for (int i = 0; i < nnz; i += nentries){
  //    for (int j = 0; j < nentries; j++){
  //       if (abs(space_values[i+j]) > 1e-14){
  //          if (space_values[i+j] < 0){
  //             values[i+j] = 1 - dtq * space_values[i+j];
  //          }
  //          else {
  //             values[i+j] = -dtq * space_values[i+j];
  //          }
  //       }
  //       else {
  //          values[i+j] = 0.0;
  //       }
  //    }
  // }


   hypre_Index *stencil_shape = hypre_StructStencilShape(stencil); 
   int diag_index = 0;
   for (int i = 0; i < hypre_StructStencilSize(stencil); i++){
      if (hypre_IndexD(stencil_shape[i], 0) == 0 && hypre_IndexD(stencil_shape[i], 1) == 0){
         diag_index = i;
      }
   }
 
   for (int i = 0; i < nnz; i += nentries){
      for (int j = 0; j < nentries; j++){
         if (j == diag_index){
            values[i+j] = 1 - dtq * space_values[i+j];
         }
         else {
            values[i+j] = -dtq * space_values[i+j];
         }
      }
   }

   //hypre_StructMatrixPrint("struct_matrix.txt", A, 1);

   HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, nentries,
                                  stencil_indices, values);
 
   HYPRE_StructMatrixAssemble(A);


   if (solver_type == SOLVER_JACOBI){
      HYPRE_StructJacobiCreate(comm, &solver);
      HYPRE_StructJacobiSetMaxIter(solver, max_iter);
      HYPRE_StructJacobiSetTol(solver, tol);
      HYPRE_StructJacobiSetup(solver, A, b, x);
      hypre_PointRelaxSetWeight(solver, 0.0);
      HYPRE_StructJacobiSolve(solver, A, b, x);
      HYPRE_StructJacobiGetNumIterations(solver, &num_iterations);
      HYPRE_StructJacobiGetFinalRelativeResidualNorm(solver, &final_res_norm);
      //printf("%d %e\n", num_iterations, final_res_norm);
      HYPRE_StructJacobiDestroy(solver);
   }
   else if (solver_type == SOLVER_JACOBI_CG){
      pcg_print_level = 0;
      HYPRE_StructPCGCreate(comm, &solver);
      HYPRE_PCGSetMaxIter((HYPRE_Solver)solver, max_iter);
      HYPRE_PCGSetTol((HYPRE_Solver)solver, tol);
      HYPRE_PCGSetTwoNorm((HYPRE_Solver)solver, 1);
      HYPRE_PCGSetRelChange((HYPRE_Solver)solver, 0);
      HYPRE_PCGSetPrintLevel((HYPRE_Solver)solver, pcg_print_level);
      HYPRE_StructJacobiCreate(comm, &precond);
      HYPRE_StructJacobiSetMaxIter(precond, 2);
      HYPRE_StructJacobiSetTol(precond, 0.0);
      HYPRE_StructJacobiSetZeroGuess(precond);
      HYPRE_PCGSetPrecond((HYPRE_Solver) solver,
                          (HYPRE_PtrToSolverFcn) HYPRE_StructJacobiSolve,
                          (HYPRE_PtrToSolverFcn) HYPRE_StructJacobiSetup,
                          (HYPRE_Solver) precond);
      HYPRE_PCGSetup((HYPRE_Solver)solver,
                     (HYPRE_Matrix)A,
                     (HYPRE_Vector)b,
                     (HYPRE_Vector)x);
      HYPRE_PCGSolve((HYPRE_Solver)solver,
                     (HYPRE_Matrix)A,
                     (HYPRE_Vector)b,
                     (HYPRE_Vector)x);
      HYPRE_PCGGetNumIterations((HYPRE_Solver)solver, &num_iterations);
      HYPRE_PCGGetFinalRelativeResidualNorm((HYPRE_Solver)solver, &final_res_norm);
      //printf("%d %e\n", num_iterations, final_res_norm);
      HYPRE_StructPCGDestroy(solver);
      HYPRE_StructJacobiDestroy(precond);

   }
   else if (solver_type == SOLVER_PFMG_CG){
      pcg_print_level = 0;
      HYPRE_StructPCGCreate(comm, &solver);
      HYPRE_PCGSetMaxIter((HYPRE_Solver)solver, max_iter);
      HYPRE_PCGSetTol((HYPRE_Solver)solver, tol);
      HYPRE_PCGSetTwoNorm((HYPRE_Solver)solver, 1);
      HYPRE_PCGSetRelChange((HYPRE_Solver)solver, 0);
      HYPRE_PCGSetPrintLevel((HYPRE_Solver)solver, pcg_print_level);
      SetupStructPFMGSolver(&precond);
      HYPRE_PCGSetPrecond((HYPRE_Solver)solver,
                          (HYPRE_PtrToSolverFcn)HYPRE_StructPFMGSolve,
                          (HYPRE_PtrToSolverFcn)HYPRE_StructPFMGSetup,
                          (HYPRE_Solver)precond);
      HYPRE_PCGSetup((HYPRE_Solver)solver,
                     (HYPRE_Matrix)A,
                     (HYPRE_Vector)b,
                     (HYPRE_Vector)x);
      HYPRE_PCGSolve((HYPRE_Solver)solver,
                     (HYPRE_Matrix)A,
                     (HYPRE_Vector)b,
                     (HYPRE_Vector)x);
      HYPRE_PCGGetNumIterations((HYPRE_Solver)solver, &num_iterations);
      HYPRE_PCGGetFinalRelativeResidualNorm((HYPRE_Solver)solver, &final_res_norm);
      //printf("%d %e\n", num_iterations, final_res_norm);
      HYPRE_StructPCGDestroy(solver);
      HYPRE_StructPFMGDestroy(precond);
   }
   else {
      SetupStructPFMGSolver(&solver);
      HYPRE_StructPFMGSetup(solver, A, b, x);
      HYPRE_StructPFMGSolve(solver, A, b, x);
      HYPRE_StructPFMGGetNumIterations(solver, &num_iterations);
      HYPRE_StructPFMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
      //printf("%d %e\n", num_iterations, final_res_norm);
      HYPRE_StructPFMGDestroy(solver);
   }

   HYPRE_StructVectorGetBoxValues(x, ilower, iupper, values);

   for (int i = 0; i < nrows; i++){
      (*y)[i] = values[i];
      (*f)[i] = ((*y)[i] - rhs[i]) / dtq;
   }
   
   HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, nentries,
                                  stencil_indices, space_values);

   free(values);
   free(space_values);
}


void HypreSolver::SetupStructPFMGSolver(HYPRE_StructSolver *pfmg_solver)
{
   HYPRE_StructPFMGCreate(comm, pfmg_solver);
   HYPRE_StructPFMGSetMaxIter(*pfmg_solver, max_iter);
   HYPRE_StructPFMGSetTol(*pfmg_solver, tol);
   HYPRE_StructPFMGSetRelChange(*pfmg_solver, 0);
   HYPRE_StructPFMGSetRAPType(*pfmg_solver, 0);
   HYPRE_StructPFMGSetRelaxType(*pfmg_solver, relax_type);
   HYPRE_StructPFMGSetJacobiWeight(*pfmg_solver, jacobi_weight);
   HYPRE_StructPFMGSetNumPreRelax(*pfmg_solver, n_pre);
   HYPRE_StructPFMGSetNumPostRelax(*pfmg_solver, n_post);
   HYPRE_StructPFMGSetSkipRelax(*pfmg_solver, 0);
   HYPRE_StructPFMGSetLogging(*pfmg_solver, 1);
}

void HypreSolver::Restrict(HYPRE_StructVector y_f, HYPRE_StructVector y_c, int f_level, int c_level)
{
   hypre_SemiRestrict(restrict_data_l[f_level], RT_l[f_level], y_f, y_c);
}

void HypreSolver::Prolong(HYPRE_StructVector y_f, HYPRE_StructVector y_c, int f_level, int c_level)
{
   hypre_SemiInterp(interp_data_l[f_level], P_l[f_level], y_c, y_f);
}

void HypreSolver::Cleanup(void)
{
   HYPRE_StructPFMGDestroy(pfmg_level_data);
   //HYPRE_StructGridDestroy(grid);
   //HYPRE_StructStencilDestroy(stencil);
   //HYPRE_StructMatrixDestroy(A);
}
