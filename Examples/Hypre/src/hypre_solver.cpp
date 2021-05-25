#include "hypre_solver.hpp"

extern int FComp_count;
extern double FComp_wtime;

HypreSolver::HypreSolver(MPI_Comm in_comm,
                         int in_dim,
                         int in_max_iter,
                         int in_max_levels,
                         int nx,
                         double in_cx, double in_cy,
                         int in_advection_flag)
{
   SetComm(in_comm);
   SetDim(in_dim);
   InitGrid(nx);

   n_pre  = 1;
   n_post = 1;
   tol = 0.0;//1e-9;
   max_iter = in_max_iter;
   jacobi_weight = 1.0;
   relax_type = RELAX_RBGS_NONSYMMETRIC;
   solver_type = SOLVER_PFMG;
   max_levels = in_max_levels;
   num_levels = in_max_levels;
   setup_done = 1;

   advection_flag = 0;
   if (in_advection_flag == 1){
      cx = in_cx;
      cy = in_cy;
      advection_flag = 1;
   }
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
                              int nx,
                              int level_index,
                              int spatial_coarsen_flag,
                              int implicit_flag,
                              double dtq)
{
   if (*A != NULL) return;

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

   if (implicit_flag == 1){
      UpdateImplicitMatrix(level_index, dtq);
   }

   HYPRE_StructMatrixAssemble(*A);
}


void HypreSolver::SetupMatrixAdvection(HYPRE_StructMatrix *A,
                                       int nx,
                                       int level_index,
                                       int spatial_coarsen_flag,
                                       int implicit_flag,
                                       double dtq)
{
   if (*A != NULL) return;

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
         values[i] = 0.0;
         for (int j = 1; j < nentries; j++){
            double c;
            if (offsets_2D[j][0] == 0){
               c = (double)offsets_2D[j][1] * cy;
            }
            else {
               c = (double)offsets_2D[j][0] * cx;
            }
            values[i+j] = c * 1.0 / h;
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

   //if (implicit_flag == 1){
   //   UpdateImplicitMatrix(level_index, dtq);
   //}

   HYPRE_StructMatrixAssemble(*A);
}

double *HypreSolver::UpdateImplicitMatrix(int level_index, double dtq)
{
   int nn = nrows_lev[level_index];
   int *iup = iupper_lev[level_index];
   int *ilow = ilower_lev[level_index];
   int nnz_l = nnz_lev[level_index];

   HYPRE_StructStencil stencil = hypre_StructMatrixStencil(A_imp_lev[level_index]);
   //int nentries = hypre_StructStencilSize(stencil);
   //hypre_Index *stencil_indices = hypre_StructStencilShape(stencil);

   double *values_poisson = (double *)calloc(nnz_l, sizeof(double));
   double *values = (double *)calloc(nnz_l, sizeof(double));

   hypre_Index *stencil_shape = hypre_StructStencilShape(stencil);
   HYPRE_StructMatrixGetBoxValues(A_imp_lev[level_index], ilow, iup, nentries, stencil_indices, values_poisson);
   int diag_index = 0;
   for (int i = 0; i < hypre_StructStencilSize(stencil); i++){
      if (hypre_IndexD(stencil_shape[i], 0) == 0 && hypre_IndexD(stencil_shape[i], 1) == 0){
         diag_index = i;
      }
   }
   for (int i = 0; i < nnz_l; i += nentries){
      for (int j = 0; j < nentries; j++){
         if (j == diag_index){
            values[i+j] = 1 - dtq * values_poisson[i+j];
         }
         else {
            values[i+j] = -dtq * values_poisson[i+j];
         }
      }
   }
   HYPRE_StructMatrixSetBoxValues(A_imp_lev[level_index], ilow, iup, nentries, stencil_indices, values);

   free(values);
   return values_poisson;
}

HYPRE_StructVector HypreSolver::SetupVectorLevel(HYPRE_StructMatrix A, int level_index)
{
   int nvalues = nrows_lev[level_index];
   int *iup = iupper_lev[level_index];
   int *ilow = ilower_lev[level_index];

   double *values = (double*) calloc(nvalues, sizeof(double));
   HYPRE_StructVector v;
   HYPRE_StructVectorCreate(comm, A->grid, &v);
   HYPRE_StructVectorInitialize(v);
   HYPRE_StructVectorSetBoxValues(v, ilow, iup, values);
   HYPRE_StructVectorAssemble(v);
   free(values);

   return v;
}

HYPRE_StructVector HypreSolver::SetupVector(void)
{
   int nvalues = nrows;
   int *iup = iupper;
   int *ilow = ilower;

   double *values = (double*) calloc(nvalues, sizeof(double));
   HYPRE_StructVector v;
   HYPRE_StructVectorCreate(comm, grid, &v);
   HYPRE_StructVectorInitialize(v);
   HYPRE_StructVectorSetBoxValues(v, ilow, iup, values);
   HYPRE_StructVectorAssemble(v);
   free(values);

   return v;
}

void HypreSolver::FEval(double *y, double t, int level_index, double **f, int piece)
{
   int nn = nrows_lev[level_index];
   if (piece != 2){
      for (int i = 0; i < nn; i++){
         (*f)[i] = 0;
      }
      return;
   }
   HYPRE_StructMatrix AA;
   //if (piece == 1){
   //   AA = A_exp_adv_lev[level_index];
   //}
   //else {
      AA = A_exp_lev[level_index];
   //}

   HYPRE_StructVector bb = b_lev[level_index];
   HYPRE_StructVector xx = x_lev[level_index];
   int *iup = iupper_lev[level_index];
   int *ilow = ilower_lev[level_index];

   double *values = (double *)calloc(nn, sizeof(double));

   int num_iterations;
   double final_res_norm;

   for (int i = 0; i < nn; i++){
      values[i] = 0.0;
   }
   HYPRE_StructVectorSetBoxValues(bb, ilow, iup, values);

   for (int i = 0; i < nn; i++){
      values[i] = y[i];
   }
   HYPRE_StructVectorSetBoxValues(xx, ilow, iup, values);

   double alpha = 1.0, beta = 0.0;
   HYPRE_StructMatrixMatvec(alpha, AA, xx, beta, bb);
   
   HYPRE_StructVectorGetBoxValues(bb, ilow, iup, values);
   for (int i = 0; i < nn; i++){
      (*f)[i] = values[i];
   }

   free(values);
}

void HypreSolver::FComp(double **y, double t, double dtq, double *rhs, int level_index, double **f)
{
   HYPRE_StructMatrix AA = A_imp_lev[level_index];
   HYPRE_StructVector bb = b_lev[level_index];
   HYPRE_StructVector xx = x_lev[level_index];
   int nn = nrows_lev[level_index];
   int nnz_l = nnz_lev[level_index];
   int *iup = iupper_lev[level_index];
   int *ilow = ilower_lev[level_index];

   int num_iterations;
   double final_res_norm;
   int setup_flag = 0;
   double *values_poisson;
   double *values = (double *)calloc(nnz_l, sizeof(double));

   for (int i = 0; i < nn; i++){
      values[i] = rhs[i];
   }
   HYPRE_StructVectorSetBoxValues(bb, ilow, iup, values);

   for (int i = 0; i < nn; i++){
      values[i] = 0.0;
   }
   HYPRE_StructVectorSetBoxValues(xx, ilow, iup, values);

   HYPRE_StructVectorAssemble(bb);
   HYPRE_StructVectorAssemble(xx);
  
   HYPRE_StructSolver hypre_solver = solver_imp_lev[level_index];
 
   if (hypre_solver == NULL){
      setup_flag = 1;
      values_poisson = UpdateImplicitMatrix(level_index, dtq);
      SetupStructSolver(level_index);
   }

   hypre_solver = solver_imp_lev[level_index];
   HYPRE_StructSolver hypre_precond = precond_imp_lev[level_index];

   if (solver_type == SOLVER_JACOBI){
      HYPRE_StructJacobiSolve(hypre_solver, A_imp_lev[level_index], bb, xx);

      //HYPRE_StructJacobiGetNumIterations(solver, &num_iterations);
      //HYPRE_StructJacobiGetFinalRelativeResidualNorm(solver, &final_res_norm);
      //printf("%d %e\n", num_iterations, final_res_norm);

      HYPRE_StructJacobiDestroy(hypre_solver);
   }
   else if (solver_type == SOLVER_JACOBI_CG){
      HYPRE_PCGSolve((HYPRE_Solver)(solver_imp_lev[level_index]),
                     (HYPRE_Matrix)(A_imp_lev[level_index]),
                     (HYPRE_Vector)bb,
                     (HYPRE_Vector)xx);

      //HYPRE_PCGGetNumIterations((HYPRE_Solver)solver, &num_iterations);
      //HYPRE_PCGGetFinalRelativeResidualNorm((HYPRE_Solver)solver, &final_res_norm);
      //printf("%d %e\n", num_iterations, final_res_norm);

      //HYPRE_StructPCGDestroy(solver);
      //HYPRE_StructJacobiDestroy(precond);

   }
   else if (solver_type == SOLVER_PFMG_CG){
      HYPRE_PCGSolve((HYPRE_Solver)(solver_imp_lev[level_index]),
                     (HYPRE_Matrix)(A_imp_lev[level_index]),
                     (HYPRE_Vector)bb,
                     (HYPRE_Vector)xx);

      //HYPRE_PCGGetNumIterations((HYPRE_Solver)solver, &num_iterations);
      //HYPRE_PCGGetFinalRelativeResidualNorm((HYPRE_Solver)solver, &final_res_norm);
      //printf("%d %e\n", num_iterations, final_res_norm);

      //HYPRE_StructPCGDestroy(solver);
      //HYPRE_StructPFMGDestroy(precond);
   }
   else {
      double wtime_start = MPI_Wtime();
      HYPRE_StructPFMGSolve(solver_imp_lev[level_index], A_imp_lev[level_index], bb, xx);
      FComp_wtime += MPI_Wtime() - wtime_start;

      //HYPRE_StructPFMGGetNumIterations(hypre_solver, &num_iterations);
      //HYPRE_StructPFMGGetFinalRelativeResidualNorm(hypre_solver, &final_res_norm);

      //printf("%d %e\n", num_iterations, final_res_norm);
   }

   HYPRE_StructVectorGetBoxValues(xx, ilow, iup, values);

   for (int i = 0; i < nn; i++){
      (*y)[i] = values[i];
      (*f)[i] = ((*y)[i] - rhs[i]) / dtq;
   }

   free(values);
   if (setup_flag == 1){
      HYPRE_StructMatrixSetBoxValues(A_imp_lev[level_index], ilow, iup, nentries, stencil_indices, values_poisson);
      free(values_poisson);
      CleanupStructSolver(&(solver_imp_lev[level_index]), &(precond_imp_lev[level_index]));
      solver_imp_lev[level_index] = NULL;
      precond_imp_lev[level_index] = NULL;
   }
}

void HypreSolver::SetupLevels(int spatial_coarsen_flag, int nx)
{
   A_imp_lev = (HYPRE_StructMatrix *)malloc(num_levels * sizeof(HYPRE_StructMatrix));
   A_exp_lev = (HYPRE_StructMatrix *)malloc(num_levels * sizeof(HYPRE_StructMatrix));
   A_exp_adv_lev = (HYPRE_StructMatrix *)malloc(num_levels * sizeof(HYPRE_StructMatrix));
   solver_imp_lev = (HYPRE_StructSolver *)malloc(num_levels * sizeof(HYPRE_StructSolver));
   precond_imp_lev = (HYPRE_StructSolver *)malloc(num_levels * sizeof(HYPRE_StructSolver));
   x_lev = (HYPRE_StructVector *)malloc(num_levels * sizeof(HYPRE_StructVector));
   b_lev = (HYPRE_StructVector *)malloc(num_levels * sizeof(HYPRE_StructVector));
   nnz_lev = (int *)malloc(num_levels * sizeof(int));
   nrows_lev = (int *)malloc(num_levels * sizeof(int));
   ilower_lev = (int **)malloc(num_levels * sizeof(int *));
   iupper_lev = (int **)malloc(num_levels * sizeof(int *));

   if (spatial_coarsen_flag == 1){
      spatial_coarsen = 1;
      SetupStructPFMGSolver(&(spatial_pfmg_data));
      HYPRE_StructPFMGSetMaxLevels(spatial_pfmg_data, 20);
      HYPRE_StructPFMGSetRAPType(spatial_pfmg_data, 0);
      HYPRE_StructPFMGSetup(spatial_pfmg_data, A_imp, b, x);
      pfmg_data = (hypre_PFMGData *)(spatial_pfmg_data);
      spatial_num_levels = pfmg_data->num_levels;
      if (spatial_num_levels < num_levels){
         printf("ERROR: Cannot further coarsen in space.  Please request fewer time levels (reduce 'nlevels' variable).\n");
         MPI_Finalize();
         exit(0);
      }
      RT_lev = pfmg_data->RT_l;
      P_lev = pfmg_data->P_l;
      x_lev = pfmg_data->x_l;
      b_lev = pfmg_data->b_l;
      restrict_data_lev = pfmg_data->restrict_data_l;
      interp_data_lev = pfmg_data->interp_data_l;

      for (int level = 0; level < num_levels; level++){
         A_imp_lev[level] = NULL;
         A_exp_lev[level] = NULL;
         A_exp_adv_lev[level] = NULL;
         solver_imp_lev[level] = NULL;
         precond_imp_lev[level] = NULL;

         int pfmg_level;
         if (level < spatial_num_levels){
            pfmg_level = level;
         }
         else {
            pfmg_level = spatial_num_levels-1;
         }

         A_imp_lev[level] = pfmg_data->A_l[pfmg_level];
         A_exp_lev[level] = pfmg_data->A_l[pfmg_level];
         //A_exp_adv_lev[level] = pfmg_data->A_l[pfmg_level];

         nrows_lev[level] = A_imp_lev[pfmg_level]->grid->local_size;
         nnz_lev[level] = nrows_lev[level] * A_imp_lev[pfmg_level]->stencil->size;
         if (nrows_lev[level] > 0){
            ilower_lev[level] = A_imp_lev[pfmg_level]->grid->boxes->boxes[0].imin;
            iupper_lev[level] = A_imp_lev[pfmg_level]->grid->boxes->boxes[0].imax;
         }
         else {
            ilower_lev[level][0] = 0;
            ilower_lev[level][1] = 0;
            iupper_lev[level][0] = -1;
            iupper_lev[level][1] = -1;
         }
         x_lev[level] = SetupVectorLevel(A_imp_lev[level], level);
         b_lev[level] = SetupVectorLevel(A_imp_lev[level], level);
      }
   }
   else {
      spatial_coarsen = 0;
      for (int level = 0; level < num_levels; level++){
         A_imp_lev[level] = NULL;
         A_exp_lev[level] = NULL;
         A_exp_adv_lev[level] = NULL;
         solver_imp_lev[level] = NULL;
         precond_imp_lev[level] = NULL;

         SetupMatrix(&(A_imp_lev[level]), nx, level, 0, 0, 0.0);
         SetupMatrix(&(A_exp_lev[level]), nx, level, 0, 0, 0.0);
         //SetupMatrixAdvection(&(A_exp_adv_lev[level]), nx, level, 0, 0, 0.0);

         nrows_lev[level] = nrows;
         nnz_lev[level] = nnz;
         ilower_lev[level] = ilower;
         iupper_lev[level] = iupper;
         x_lev[level] = SetupVectorLevel(A_imp_lev[level], level);
         b_lev[level] = SetupVectorLevel(A_imp_lev[level], level);
      }
   }
}

void HypreSolver::SetupStructSolver(int level_index)
{
   HYPRE_StructMatrix AA = A_imp_lev[level_index];
   HYPRE_StructVector bb = b_lev[level_index];
   HYPRE_StructVector xx = x_lev[level_index];   
   HYPRE_StructSolver hypre_solver = NULL;
   HYPRE_StructSolver hypre_precond = NULL;

   if (solver_type == SOLVER_JACOBI){
      HYPRE_StructJacobiCreate(comm, &hypre_solver);
      HYPRE_StructJacobiSetMaxIter(hypre_solver, max_iter);
      HYPRE_StructJacobiSetTol(hypre_solver, tol);
      HYPRE_StructJacobiSetup(hypre_solver, AA, bb, xx);
   }
   else if (solver_type == SOLVER_JACOBI_CG){
      pcg_print_level = 0;
      HYPRE_StructPCGCreate(comm, &hypre_solver);
      HYPRE_PCGSetMaxIter((HYPRE_Solver)hypre_solver, max_iter);
      HYPRE_PCGSetTol((HYPRE_Solver)hypre_solver, tol);
      HYPRE_PCGSetTwoNorm((HYPRE_Solver)hypre_solver, 1);
      HYPRE_PCGSetRelChange((HYPRE_Solver)hypre_solver, 0);
      HYPRE_PCGSetPrintLevel((HYPRE_Solver)hypre_solver, pcg_print_level);
      HYPRE_StructJacobiCreate(comm, &hypre_precond);
      HYPRE_StructJacobiSetMaxIter(hypre_precond, 2);
      HYPRE_StructJacobiSetTol(hypre_precond, 0.0);
      HYPRE_StructJacobiSetZeroGuess(hypre_precond);
      HYPRE_PCGSetPrecond((HYPRE_Solver)hypre_solver,
                          (HYPRE_PtrToSolverFcn) HYPRE_StructJacobiSolve,
                          (HYPRE_PtrToSolverFcn) HYPRE_StructJacobiSetup,
                          (HYPRE_Solver)hypre_precond);
      HYPRE_PCGSetup((HYPRE_Solver)hypre_solver,
                     (HYPRE_Matrix)AA,
                     (HYPRE_Vector)bb,
                     (HYPRE_Vector)xx);

   }
   else if (solver_type == SOLVER_PFMG_CG){
      pcg_print_level = 0;
      HYPRE_StructPCGCreate(comm, &hypre_solver);
      HYPRE_PCGSetMaxIter((HYPRE_Solver)hypre_solver, max_iter);
      HYPRE_PCGSetTol((HYPRE_Solver)hypre_solver, tol);
      HYPRE_PCGSetTwoNorm((HYPRE_Solver)hypre_solver, 1);
      HYPRE_PCGSetRelChange((HYPRE_Solver)hypre_solver, 0);
      HYPRE_PCGSetPrintLevel((HYPRE_Solver)hypre_solver, pcg_print_level);
      SetupStructPFMGSolver(&hypre_precond);
      HYPRE_PCGSetPrecond((HYPRE_Solver)hypre_solver,
                          (HYPRE_PtrToSolverFcn)HYPRE_StructPFMGSolve,
                          (HYPRE_PtrToSolverFcn)HYPRE_StructPFMGSetup,
                          (HYPRE_Solver)hypre_precond);
      HYPRE_PCGSetup((HYPRE_Solver)hypre_solver,
                     (HYPRE_Matrix)AA,
                     (HYPRE_Vector)bb,
                     (HYPRE_Vector)xx);
   }
   else {
      SetupStructPFMGSolver(&hypre_solver);
      HYPRE_StructPFMGSetup(hypre_solver, AA, bb, xx);
   }

   solver_imp_lev[level_index] = hypre_solver;
   precond_imp_lev[level_index] = hypre_precond;
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
   HYPRE_StructPFMGSetSkipRelax(*pfmg_solver, 1);
   HYPRE_StructPFMGSetLogging(*pfmg_solver, 1);
   HYPRE_StructPFMGSetPrintLevel(*pfmg_solver, 1);
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
   if (solver_imp != NULL) CleanupStructSolver(&solver_imp, &precond_imp);
   if (A_imp != NULL) HYPRE_StructMatrixDestroy(A_imp);
   if (A_exp != NULL) HYPRE_StructMatrixDestroy(A_exp);
   if (x != NULL) HYPRE_StructVectorDestroy(x);
   if (b != NULL) HYPRE_StructVectorDestroy(b);

 
   //if (spatial_pfmg_data != NULL) HYPRE_StructPFMGDestroy(spatial_pfmg_data);
   if (solver_imp_lev != NULL) {
      for (int level_index = 0; level_index < num_levels; level_index++){
         if (solver_imp_lev[level_index] != NULL)
            CleanupStructSolver(&(solver_imp_lev[level_index]), &(precond_imp_lev[level_index]));
      }
      free(solver_imp_lev);
      free(precond_imp_lev);
   } 
   if (A_imp_lev != NULL) {
      for (int level_index = 0; level_index < num_levels; level_index++){
         if (A_imp_lev[level_index] != NULL) HYPRE_StructMatrixDestroy(A_imp_lev[level_index]);
      }
      free(A_imp_lev);
   }
   if (A_exp_lev != NULL) {
      if (spatial_coarsen == 0){
         for (int level_index = 0; level_index < num_levels; level_index++){
            if (A_exp_lev[level_index] != NULL) HYPRE_StructMatrixDestroy(A_exp_lev[level_index]);
         }
      }
      free(A_exp_lev);
   }
   if (b_lev != NULL) {
      for (int level_index = 0; level_index < num_levels; level_index++){
         if (b_lev[level_index] != NULL) HYPRE_StructVectorDestroy(b_lev[level_index]);
      }
      free(b_lev);
   }
   if (x_lev != NULL) {
      for (int level_index = 0; level_index < num_levels; level_index++){
         if (x_lev[level_index] != NULL) HYPRE_StructVectorDestroy(x_lev[level_index]);
      }
      free(x_lev);
   }
   if (nnz_lev != NULL) free(nnz_lev);
   if (nrows_lev != NULL) free(nrows_lev);
   if (ilower_lev != NULL){
      for (int level_index = 0; level_index < num_levels; level_index++){
         if (ilower_lev[level_index] != NULL) free(ilower_lev[level_index]);
      }
      free(ilower_lev);
   }
   if (iupper_lev != NULL){
      for (int level_index = 0; level_index < num_levels; level_index++){
         if (iupper_lev[level_index] != NULL) free(iupper_lev[level_index]);
      }
      free(iupper_lev);
   }
}
