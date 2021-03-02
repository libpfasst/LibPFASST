#include "hypre_solver.hpp"

HypreSolver::HypreSolver(MPI_Comm in_comm, int in_dim, int in_max_iter, int in_max_levels, int num_grid_points)
{
   //if (setup_done == 1) return;
   SetComm(in_comm);
   SetDim(in_dim);
   InitGrid(num_grid_points);

   n_pre  = 1;
   n_post = 1;
   tol = 1e-9;
   max_iter = in_max_iter;
   jacobi_weight = 1.0;
   relax_type = RELAX_RBGS;
   solver_type = SOLVER_PFMG;
   max_levels = in_max_levels;
   num_levels = in_max_levels;
   setup_done = 1;
}

HypreSolver::~HypreSolver(void)
{
   Cleanup();
}

int HypreSolver::GetNumRowsLevel(int level_index)
{
   return nrows_lev[level_index];
}

int HypreSolver::GetNumLevels(void)
{
   return num_levels;
}

void HypreSolver::SetupMatrix(HYPRE_StructMatrix *A,
                              int num_grid_points,
                              int level_index,
                              int spacial_coarsen_flag,
                              int implicit_flag,
                              double dtq)
{
   if ((level_index > 0) && (spacial_coarsen_flag == 1)) return;

   if (*A == NULL){
      double *values;

      /* Set up a Struct Matrix */
      {
         int nvalues = nnz;

         /* Create an empty matrix object */
         HYPRE_StructMatrixCreate(comm, grid, stencil, A);

         /* Indicate that the matrix coefficients are ready to be set */
         HYPRE_StructMatrixInitialize(*A);

         values = (double *)calloc(nvalues, sizeof(double));

         for (int j = 0; j < nentries; j++){
            stencil_indices[j] = j;
         }

         /* Set the standard stencil at each grid point,
            we will fix the boundaries later */
         for (int i = 0; i < nvalues; i += nentries){
            values[i] = -4.0 / h2;
            for (int j = 1; j < nentries; j++){
               values[i+j] = 1.0 / h2;
            }
         }

         HYPRE_StructMatrixSetBoxValues(*A, ilower, iupper, nentries,
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
         /* Recall: pi and pj describe position in the processor grid */
         if (pj == 0){
            /* Bottom row of grid points */
            bc_ilower[0] = pi*n;
            bc_ilower[1] = pj*n;

            bc_iupper[0] = bc_ilower[0] + n-1;
            bc_iupper[1] = bc_ilower[1];

            bc_stencil_indices[0] = 3;

            HYPRE_StructMatrixSetBoxValues(*A, bc_ilower, bc_iupper, bc_nentries,
                                           bc_stencil_indices, values);
         }

         if (pj == N-1){
            /* upper row of grid points */
            bc_ilower[0] = pi*n;
            bc_ilower[1] = pj*n + n-1;

            bc_iupper[0] = bc_ilower[0] + n-1;
            bc_iupper[1] = bc_ilower[1];

            bc_stencil_indices[0] = 4;

            HYPRE_StructMatrixSetBoxValues(*A, bc_ilower, bc_iupper, bc_nentries,
                                           bc_stencil_indices, values);
         }

         if (pi == 0){
            /* Left row of grid points */
            bc_ilower[0] = pi*n;
            bc_ilower[1] = pj*n;

            bc_iupper[0] = bc_ilower[0];
            bc_iupper[1] = bc_ilower[1] + n-1;

            bc_stencil_indices[0] = 1;

            HYPRE_StructMatrixSetBoxValues(*A, bc_ilower, bc_iupper, bc_nentries,
                                           bc_stencil_indices, values);
         }

         if (pi == N-1){
            /* Right row of grid points */
            bc_ilower[0] = pi*n + n-1;
            bc_ilower[1] = pj*n;

            bc_iupper[0] = bc_ilower[0];
            bc_iupper[1] = bc_ilower[1] + n-1;

            bc_stencil_indices[0] = 2;

            HYPRE_StructMatrixSetBoxValues(*A, bc_ilower, bc_iupper, bc_nentries,
                                           bc_stencil_indices, values);
         }
         free(values);
      }
   }

   if (implicit_flag == 1){
      UpdateImplicitMatrix(A, dtq);
   }

   HYPRE_StructMatrixAssemble(*A);
}

double* HypreSolver::UpdateImplicitMatrix(HYPRE_StructMatrix *A, double dtq)
{
   double *values_poisson = (double *)calloc(nnz, sizeof(double));
   double *values = (double *)calloc(nnz, sizeof(double));

   HYPRE_StructMatrixGetBoxValues(*A, ilower, iupper, nentries,
                                  stencil_indices, values_poisson);
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
            values[i+j] = 1 - dtq * values_poisson[i+j];
         }
         else {
            values[i+j] = -dtq * values_poisson[i+j];
         }
      }
   }
   HYPRE_StructMatrixSetBoxValues(*A, ilower, iupper, nentries,
                                  stencil_indices, values);

   free(values);
   return values_poisson;
}

HYPRE_StructVector HypreSolver::SetupVector(void)
{
   int nvalues = nrows;
   double *values = (double*) calloc(nvalues, sizeof(double));
   HYPRE_StructVector v;
   HYPRE_StructVectorCreate(comm, grid, &v);
   HYPRE_StructVectorInitialize(v);
   HYPRE_StructVectorSetBoxValues(v, ilower, iupper, values);
   HYPRE_StructVectorAssemble(v);
   free(values);
   return v;
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

   double alpha = 1.0, beta = 0.0;
   HYPRE_StructMatrixMatvec(alpha, A_exp, x, beta, b);
   

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
   int setup_flag = 0;
   double *values_poisson;
   double *values = (double *)calloc(nnz, sizeof(double));

   for (int i = 0; i < nrows; i++){
      values[i] = rhs[i];
   }
   HYPRE_StructVectorSetBoxValues(b, ilower, iupper, values);

   for (int i = 0; i < nrows; i++){
      values[i] = 0.0;
   }
   HYPRE_StructVectorSetBoxValues(x, ilower, iupper, values);

   //HYPRE_StructVectorAssemble(b);
   //HYPRE_StructVectorAssemble(x);
   
   if (solver_imp == NULL){
      setup_flag = 1;
      values_poisson = UpdateImplicitMatrix(&A_imp, dtq);
      SetupStructSolver(&A_imp, &solver_imp, &precond_imp);
   }

   if (solver_type == SOLVER_JACOBI){
      HYPRE_StructJacobiSolve(solver_imp, A_imp, b, x);

      //HYPRE_StructJacobiGetNumIterations(solver, &num_iterations);
      //HYPRE_StructJacobiGetFinalRelativeResidualNorm(solver, &final_res_norm);
      //printf("%d %e\n", num_iterations, final_res_norm);

      HYPRE_StructJacobiDestroy(solver_imp);
   }
   else if (solver_type == SOLVER_JACOBI_CG){
      HYPRE_PCGSolve((HYPRE_Solver)solver_imp,
                     (HYPRE_Matrix)A_imp,
                     (HYPRE_Vector)b,
                     (HYPRE_Vector)x);

      //HYPRE_PCGGetNumIterations((HYPRE_Solver)solver, &num_iterations);
      //HYPRE_PCGGetFinalRelativeResidualNorm((HYPRE_Solver)solver, &final_res_norm);
      //printf("%d %e\n", num_iterations, final_res_norm);

      //HYPRE_StructPCGDestroy(solver);
      //HYPRE_StructJacobiDestroy(precond);

   }
   else if (solver_type == SOLVER_PFMG_CG){
      HYPRE_PCGSolve((HYPRE_Solver)solver_imp,
                     (HYPRE_Matrix)A_imp,
                     (HYPRE_Vector)b,
                     (HYPRE_Vector)x);

      //HYPRE_PCGGetNumIterations((HYPRE_Solver)solver, &num_iterations);
      //HYPRE_PCGGetFinalRelativeResidualNorm((HYPRE_Solver)solver, &final_res_norm);
      //printf("%d %e\n", num_iterations, final_res_norm);

      //HYPRE_StructPCGDestroy(solver);
      //HYPRE_StructPFMGDestroy(precond);
   }
   else {
      HYPRE_StructPFMGSolve(solver_imp, A_imp, b, x);

      //HYPRE_StructPFMGGetNumIterations(solver_imp, &num_iterations);
      //HYPRE_StructPFMGGetFinalRelativeResidualNorm(solver, &final_res_norm);

      //printf("%d %e\n", num_iterations, final_res_norm);
      //HYPRE_StructPFMGDestroy(solver);
   }

   HYPRE_StructVectorGetBoxValues(x, ilower, iupper, values);

   for (int i = 0; i < nrows; i++){
      (*y)[i] = values[i];
      (*f)[i] = ((*y)[i] - rhs[i]) / dtq;
   }
   

   free(values);
   if (setup_flag == 1){
      HYPRE_StructMatrixSetBoxValues(A_imp, ilower, iupper, nentries,
                                     stencil_indices, values_poisson);
      free(values_poisson);
      CleanupStructSolver(&solver_imp, &precond_imp);
      solver_imp = NULL;
   }
}

void HypreSolver::SetupSpacialCoarsen(void)
{
   //SetupStructPFMGSolver(&(pfmg_level_data));
   //HYPRE_StructPFMGSetMaxLevels(pfmg_level_data, max_levels);
   //HYPRE_StructPFMGSetup(pfmg_level_data, A, b, x);
   //pfmg_data = (hypre_PFMGData *)(pfmg_level_data);
   //A_lev = pfmg_data->A_l;
   //RT_lev = pfmg_data->RT_l;
   //P_lev = pfmg_data->P_l;
   //x_lev = pfmg_data->x_l;
   //b_lev = pfmg_data->b_l;
   //restrict_data_lev = pfmg_data->restrict_data_l;
   //interp_data_lev = pfmg_data->interp_data_l;
   //num_levels = pfmg_data->num_levels;

   //nrows_lev = (int *)malloc(num_levels * sizeof(int));
   //ilower_lev = (int **)malloc(num_levels * sizeof(int *));
   //iupper_lev = (int **)malloc(num_levels * sizeof(int *));
   //for (int level = 0; level < num_levels; level++){
   //   ilower_lev[level] = (int *)malloc(2 * sizeof(int));
   //   iupper_lev[level] = (int *)malloc(2 * sizeof(int));
   //   nrows_lev[level] = A_lev[level]->grid->local_size;
   //   if (nrows_l[level] > 0){
   //      ilower_lev[level] = grid->boxes->boxes[0].imin;
   //      iupper_lev[level] = grid->boxes->boxes[0].imax;
   //   }
   //   else {
   //      ilower_lev[level][0] = 0;
   //      ilower_lev[level][1] = 0;
   //      iupper_lev[level][0] = -1;
   //      iupper_lev[level][1] = -1;
   //   }
   //}
}

void HypreSolver::SetupStructSolver(HYPRE_StructMatrix *A, HYPRE_StructSolver *solver, HYPRE_StructSolver *precond) 
{
   if (solver_type == SOLVER_JACOBI){
      HYPRE_StructJacobiCreate(comm, solver);
      HYPRE_StructJacobiSetMaxIter(*solver, max_iter);
      HYPRE_StructJacobiSetTol(*solver, tol);
      HYPRE_StructJacobiSetup(*solver, *A, b, x);
   }
   else if (solver_type == SOLVER_JACOBI_CG){
      pcg_print_level = 0;
      HYPRE_StructPCGCreate(comm, solver);
      HYPRE_PCGSetMaxIter((HYPRE_Solver)(*solver), max_iter);
      HYPRE_PCGSetTol((HYPRE_Solver)(*solver), tol);
      HYPRE_PCGSetTwoNorm((HYPRE_Solver)(*solver), 1);
      HYPRE_PCGSetRelChange((HYPRE_Solver)(*solver), 0);
      HYPRE_PCGSetPrintLevel((HYPRE_Solver)(*solver), pcg_print_level);
      HYPRE_StructJacobiCreate(comm, precond);
      HYPRE_StructJacobiSetMaxIter(*precond, 2);
      HYPRE_StructJacobiSetTol(*precond, 0.0);
      HYPRE_StructJacobiSetZeroGuess(*precond);
      HYPRE_PCGSetPrecond((HYPRE_Solver) (*solver),
                          (HYPRE_PtrToSolverFcn) HYPRE_StructJacobiSolve,
                          (HYPRE_PtrToSolverFcn) HYPRE_StructJacobiSetup,
                          (HYPRE_Solver) (*precond));
      HYPRE_PCGSetup((HYPRE_Solver)(*solver),
                     (HYPRE_Matrix)(*A),
                     (HYPRE_Vector)b,
                     (HYPRE_Vector)x);

   }
   else if (solver_type == SOLVER_PFMG_CG){
      pcg_print_level = 0;
      HYPRE_StructPCGCreate(comm, solver);
      HYPRE_PCGSetMaxIter((HYPRE_Solver)(*solver), max_iter);
      HYPRE_PCGSetTol((HYPRE_Solver)(*solver), tol);
      HYPRE_PCGSetTwoNorm((HYPRE_Solver)(*solver), 1);
      HYPRE_PCGSetRelChange((HYPRE_Solver)(*solver), 0);
      HYPRE_PCGSetPrintLevel((HYPRE_Solver)(*solver), pcg_print_level);
      SetupStructPFMGSolver(precond);
      HYPRE_PCGSetPrecond((HYPRE_Solver)(*solver),
                          (HYPRE_PtrToSolverFcn)HYPRE_StructPFMGSolve,
                          (HYPRE_PtrToSolverFcn)HYPRE_StructPFMGSetup,
                          (HYPRE_Solver)(*precond));
      HYPRE_PCGSetup((HYPRE_Solver)(*solver),
                     (HYPRE_Matrix)(*A),
                     (HYPRE_Vector)b,
                     (HYPRE_Vector)x);
   }
   else {
      SetupStructPFMGSolver(solver);
      HYPRE_StructPFMGSetup(*solver, *A, b, x);
   }
}

void HypreSolver::CleanupStructSolver(HYPRE_StructSolver *solver, HYPRE_StructSolver *precond)
{
   if (solver_type == SOLVER_JACOBI){
      HYPRE_StructJacobiDestroy(*solver);
   }
   else if (solver_type == SOLVER_JACOBI_CG){
      HYPRE_StructPCGDestroy(*solver);
      HYPRE_StructJacobiDestroy(*precond);
   }
   else if (solver_type == SOLVER_PFMG_CG){
      HYPRE_StructPCGDestroy(*solver);
      HYPRE_StructPFMGDestroy(*precond);
   }
   else {
      HYPRE_StructPFMGDestroy(*solver);
   }
}

void HypreSolver::SetupStructPFMGSolver(HYPRE_StructSolver *pfmg_solver)
{
   HYPRE_StructPFMGCreate(comm, pfmg_solver);
   HYPRE_StructPFMGSetMaxIter(*pfmg_solver, max_iter);
   HYPRE_StructPFMGSetTol(*pfmg_solver, tol);
   HYPRE_StructPFMGSetRelChange(*pfmg_solver, 0);
   HYPRE_StructPFMGSetRAPType(*pfmg_solver, 1);
   HYPRE_StructPFMGSetRelaxType(*pfmg_solver, relax_type);
   HYPRE_StructPFMGSetJacobiWeight(*pfmg_solver, jacobi_weight);
   HYPRE_StructPFMGSetNumPreRelax(*pfmg_solver, n_pre);
   HYPRE_StructPFMGSetNumPostRelax(*pfmg_solver, n_post);
   HYPRE_StructPFMGSetSkipRelax(*pfmg_solver, 0);
   HYPRE_StructPFMGSetLogging(*pfmg_solver, 1);
}

void HypreSolver::Restrict(HYPRE_StructVector y_f, HYPRE_StructVector y_c, int f_level, int c_level)
{
   hypre_SemiRestrict(restrict_data_lev[f_level], RT_lev[f_level], y_f, y_c);
}

void HypreSolver::Prolong(HYPRE_StructVector y_f, HYPRE_StructVector y_c, int f_level, int c_level)
{
   hypre_SemiInterp(interp_data_lev[f_level], P_lev[f_level], y_c, y_f);
}

void HypreSolver::Cleanup(void)
{
   if (grid != NULL) HYPRE_StructGridDestroy(grid);
   if (stencil != NULL) HYPRE_StructStencilDestroy(stencil);
   if (A_imp != NULL) HYPRE_StructMatrixDestroy(A_imp);
   if (A_exp != NULL) HYPRE_StructMatrixDestroy(A_exp);
   if (solver_imp != NULL) CleanupStructSolver(&solver_imp, &precond_imp);
   if (b != NULL) HYPRE_StructVectorDestroy(b);
   if (x != NULL) HYPRE_StructVectorDestroy(x);
}
