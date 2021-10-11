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
   InitGrid(nx, -1, NULL, NULL, NULL, 1);

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

      HYPRE_StructMatrixSetBoxValues(*A, ilower, iupper, nentries, stencil_indices, values);

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

         HYPRE_StructMatrixSetBoxValues(*A, bc_ilower, bc_iupper, bc_nentries, bc_stencil_indices, values);
      }

      if (pj == N-1){
         /* upper row of grid points */
         bc_ilower[0] = pi*n;
         bc_ilower[1] = pj*n + n-1;

         bc_iupper[0] = bc_ilower[0] + n-1;
         bc_iupper[1] = bc_ilower[1];

         bc_stencil_indices[0] = 4;

         HYPRE_StructMatrixSetBoxValues(*A, bc_ilower, bc_iupper, bc_nentries, bc_stencil_indices, values);
      }

      if (pi == 0){
         /* Left row of grid points */
         bc_ilower[0] = pi*n;
         bc_ilower[1] = pj*n;

         bc_iupper[0] = bc_ilower[0];
         bc_iupper[1] = bc_ilower[1] + n-1;

         bc_stencil_indices[0] = 1;

         HYPRE_StructMatrixSetBoxValues(*A, bc_ilower, bc_iupper, bc_nentries, bc_stencil_indices, values);
      }

      if (pi == N-1){
         /* Right row of grid points */
         bc_ilower[0] = pi*n + n-1;
         bc_ilower[1] = pj*n;

         bc_iupper[0] = bc_ilower[0];
         bc_iupper[1] = bc_ilower[1] + n-1;

         bc_stencil_indices[0] = 2;

         HYPRE_StructMatrixSetBoxValues(*A, bc_ilower, bc_iupper, bc_nentries, bc_stencil_indices, values);
      }
      free(values);
   }

   if (implicit_flag == 1){
      UpdateImplicitMatrix(level_index, dtq);
   }

   HYPRE_StructMatrixAssemble(*A);
}

double *HypreSolver::UpdateImplicitMatrix(int level_index, double dtq)
{
   int nrows_l = nrows_lev[level_index];
   int *iup = iupper_lev[level_index];
   int *ilow = ilower_lev[level_index];
   int nnz_l = nnz_lev[level_index];

   HYPRE_StructStencil stencil = hypre_StructMatrixStencil(A_imp_lev[level_index]);
   int stencil_size = hypre_StructStencilSize(stencil);

   double *values_poisson = (double *)calloc(nnz_l, sizeof(double));
   double *values = (double *)calloc(nnz_l, sizeof(double));
   int *stencil_ind = (int *)calloc(stencil_size, sizeof(int));
   for (int j = 0; j < stencil_size; j++) stencil_ind[j] = j;

   hypre_Index *stencil_shape = hypre_StructStencilShape(stencil);
   HYPRE_StructMatrixGetBoxValues(A_imp_lev[level_index], ilow, iup, stencil_size, stencil_ind, values_poisson);
   int diag_index = 0;
   for (int i = 0; i < stencil_size; i++){
      if (hypre_IndexD(stencil_shape[i], 0) == 0 && hypre_IndexD(stencil_shape[i], 1) == 0){
         diag_index = i;
      }
   }
   for (int i = 0; i < nnz_l; i += stencil_size){
      for (int j = 0; j < stencil_size; j++){
         if (j == diag_index){
            values[i+j] = 1 - dtq * values_poisson[i+j];
         }
         else {
            values[i+j] = -dtq * values_poisson[i+j];
         }
      }
   }
   HYPRE_StructMatrixSetBoxValues(A_imp_lev[level_index], ilow, iup, stencil_size, stencil_ind, values);

   //int level = level_index;
   //char filename[100];
   //sprintf(filename, "../../../matlab_codes/A_imp_%d", level);
   //HYPRE_StructMatrixPrint(filename, A_imp_lev[level], 0);

   free(stencil_ind);
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
   HYPRE_StructVectorCreate(comm, hypre_StructMatrixGrid(A), &v);
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
   int nrows_l = nrows_lev[level_index];
   if (piece != 2){
      for (int i = 0; i < nrows_l; i++){
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

   double *values = (double *)calloc(nrows_l, sizeof(double));

   int num_iterations;
   double final_res_norm;

   for (int i = 0; i < nrows_l; i++){
      values[i] = 0.0;
   }
   HYPRE_StructVectorSetBoxValues(bb, ilow, iup, values);

   for (int i = 0; i < nrows_l; i++){
      values[i] = y[i];
   }
   HYPRE_StructVectorSetBoxValues(xx, ilow, iup, values);

   double alpha = 1.0, beta = 0.0;
   HYPRE_StructMatrixMatvec(alpha, AA, xx, beta, bb);
   
   HYPRE_StructVectorGetBoxValues(bb, ilow, iup, values);
   for (int i = 0; i < nrows_l; i++){
      (*f)[i] = values[i];
   }

   free(values);
}

void HypreSolver::FComp(double **y, double t, double dtq, double *rhs, int level_index, double **f)
{
   HYPRE_StructMatrix AA = A_imp_lev[level_index];
   HYPRE_StructVector bb = b_lev[level_index];
   HYPRE_StructVector xx = x_lev[level_index];
   int nrows_l = nrows_lev[level_index];
   int nnz_l = nnz_lev[level_index];
   int *iup = iupper_lev[level_index];
   int *ilow = ilower_lev[level_index];

   int num_iterations;
   double final_res_norm;
   int setup_flag = 0;
   double *values_poisson;
   double *values = (double *)calloc(nnz_l, sizeof(double));

   for (int i = 0; i < nrows_l; i++){
      values[i] = rhs[i];
   }
   HYPRE_StructVectorSetBoxValues(bb, ilow, iup, values);

   for (int i = 0; i < nrows_l; i++){
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

   for (int i = 0; i < nrows_l; i++){
      (*y)[i] = values[i];
      (*f)[i] = ((*y)[i] - rhs[i]) / dtq;
   }

   free(values);
   if (setup_flag == 1){      
      HYPRE_StructStencil stencil = hypre_StructMatrixStencil(A_imp_lev[level_index]);
      int stencil_size = hypre_StructStencilSize(stencil);
      int *stencil_ind = (int *)calloc(stencil_size, sizeof(int));
      for (int j = 0; j < stencil_size; j++) stencil_ind[j] = j;
      HYPRE_StructMatrixSetBoxValues(A_imp_lev[level_index], ilow, iup, stencil_size, stencil_ind, values_poisson);
      free(stencil_ind);
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

      P_lev = (HYPRE_StructMatrix *)malloc(num_levels * sizeof(HYPRE_StructMatrix));
      R_lev = (HYPRE_StructMatrix *)malloc(num_levels * sizeof(HYPRE_StructMatrix));

      HYPRE_StructPFMGCreate(comm, &coarse_space_data);
      HYPRE_StructPFMGSetMaxLevels(coarse_space_data, 20);
      HYPRE_StructPFMGSetRelaxType(coarse_space_data, RELAX_JACOBI);
      HYPRE_StructPFMGSetJacobiWeight(coarse_space_data, .5);
      HYPRE_StructPFMGSetRAPType(coarse_space_data, 1);
      HYPRE_StructPFMGSetSkipRelax(coarse_space_data, 0);
      HYPRE_StructPFMGSetup(coarse_space_data, A_imp, b, x);
      hypre_PFMGData *pfmg_data = (hypre_PFMGData *)(coarse_space_data);
      spatial_num_levels = pfmg_data->num_levels;
      //R_lev = pfmg_data->RT_l;
      //P_lev = pfmg_data->P_l;
      restrict_data_lev = pfmg_data->restrict_data_l;
      interp_data_lev = pfmg_data->interp_data_l;

      if (spatial_num_levels < num_levels){
         printf("ERROR: Cannot further coarsen in space.  Please request fewer time levels (reduce 'nlevels' variable).\n");
         MPI_Finalize();
         exit(0);
      }

      for (int level = 0; level < num_levels; level++){
         A_imp_lev[level] = NULL;
         A_exp_lev[level] = NULL;
         A_exp_adv_lev[level] = NULL;
         solver_imp_lev[level] = NULL;
         precond_imp_lev[level] = NULL;
         x_lev[level] = NULL;
         b_lev[level] = NULL;
         P_lev[level] = NULL;
         R_lev[level] = NULL;

         int space_level;
         if (level < spatial_num_levels){
            space_level = level;
         }
         else {
            space_level = spatial_num_levels-1;
         }

         HYPRE_StructMatrix A_l = pfmg_data->A_l[space_level];
         HYPRE_StructStencil stencil = hypre_StructMatrixStencil(A_l);
         HYPRE_StructGrid grid = hypre_StructMatrixGrid(A_l);
         int stencil_size = hypre_StructStencilSize(stencil);

         nrows_lev[level] = hypre_StructGridLocalSize(grid);
         if (level > 0){
            int min_nrows; 
            MPI_Allreduce(&nrows_lev[level], &min_nrows, 1, MPI_INT, MPI_MIN, comm);
            if (min_nrows < 1){
               printf("ERROR: Cannot further coarsen in space.  Please request fewer time levels (reduce 'nlevels' variable).\n");
               MPI_Finalize();
               exit(0);
            }
         }
         
         nnz_lev[level] = nrows_lev[level] * stencil_size;
         if (nrows_lev[level] > 0){
            ilower_lev[level] = grid->boxes->boxes[0].imin;
            iupper_lev[level] = grid->boxes->boxes[0].imax;
         }
         else {
            ilower_lev[level][0] = 0;
            ilower_lev[level][1] = 0;
            iupper_lev[level][0] = -1;
            iupper_lev[level][1] = -1;
         }

         int nrows_l, nnz_l, *ilow, *iup, *stencil_ind;
         double *values;         

         nrows_l = hypre_StructGridLocalSize(grid);
         nnz_l = nrows_l * stencil_size;
         values = (double *)calloc(nnz_l, sizeof(double));
         ilow = grid->boxes->boxes[0].imin;
         iup = grid->boxes->boxes[0].imax;
         stencil_ind = (int *)calloc(stencil_size, sizeof(int));
         for (int j = 0; j < stencil_size; j++) stencil_ind[j] = j;

         HYPRE_StructMatrixCreate(comm,
                                  grid,
                                  stencil,
                                  &(A_imp_lev[level]));
         HYPRE_StructMatrixInitialize(A_imp_lev[level]);
         HYPRE_StructMatrixGetBoxValues(A_l, ilow, iup, stencil_size, stencil_ind, values);
         hypre_Box *box = hypre_BoxArrayBox(hypre_StructGridBoxes(grid), 0);
         int nx = hypre_BoxIMaxD(box, 0);
         int ny = hypre_BoxIMaxD(box, 1);
         double hx2 = pow(1.0/(double)(nx+2), 2.0);
         double hy2 = pow(1.0/(double)(ny+2), 2.0);
         hypre_Index *stencil_shape = hypre_StructStencilShape(stencil);
         int diag_index = 0;
         for (int i = 0; i < stencil_size; i++){
            if (hypre_IndexD(stencil_shape[i], 0) == 0 && hypre_IndexD(stencil_shape[i], 1) == 0){
               diag_index = i;
            }
         }
         for (int i = 0; i < nnz_l; i += stencil_size){
            for (int j = 0; j < stencil_size; j++){
               if (j == diag_index){
                  values[i+j] = -(2.0 / hx2 + 2.0 / hy2);
               }
               else {
                  if (hypre_IndexD(stencil_shape[j], 0) == 0){
                     values[i+j] = 1.0 / hy2;
                  }
                  else {
                     values[i+j] = 1.0 / hx2;
                  }
               }
            }
         }
         HYPRE_StructMatrixSetBoxValues(A_imp_lev[level], ilow, iup, stencil_size, stencil_ind, values);
         hypre_StructMatrixClearBoundary(A_imp_lev[level]);
         HYPRE_StructMatrixAssemble(A_imp_lev[level]);

         HYPRE_StructMatrixCreate(comm,
                                  grid,
                                  stencil,
                                  &(A_exp_lev[level]));
         HYPRE_StructMatrixInitialize(A_exp_lev[level]);
         HYPRE_StructMatrixSetBoxValues(A_exp_lev[level], ilow, iup, stencil_size, stencil_ind, values);
         HYPRE_StructMatrixAssemble(A_exp_lev[level]);

         free(stencil_ind);
         free(values);

         if (space_level < spatial_num_levels-1){
            HYPRE_StructMatrix P_l = pfmg_data->P_l[space_level];
            stencil = hypre_StructMatrixStencil(P_l);
            grid = hypre_StructMatrixGrid(P_l);
            stencil_size = hypre_StructStencilSize(stencil);
            nrows_l = hypre_StructGridLocalSize(grid);
            nnz_l = nrows_l * stencil_size;
            values = (double *)calloc(nnz_l, sizeof(double));
            ilow = grid->boxes->boxes[0].imin;
            iup = grid->boxes->boxes[0].imax;
            stencil_ind = (int *)calloc(stencil_size, sizeof(int));
            for (int j = 0; j < stencil_size; j++) stencil_ind[j] = j;

            HYPRE_StructMatrixCreate(comm,
                                     grid,
                                     stencil,
                                     &(P_lev[level]));
            HYPRE_StructMatrixInitialize(P_lev[level]);
            HYPRE_StructMatrixGetBoxValues(P_l, ilow, iup, stencil_size, stencil_ind, values);
            for (int i = 0; i < nnz_l; i++) values[i] = .5;
            HYPRE_StructMatrixSetBoxValues(P_lev[level], ilow, iup, stencil_size, stencil_ind, values);
            //hypre_StructMatrixClearBoundary(P_lev[level]);
            HYPRE_StructMatrixAssemble(P_lev[level]);

            free(stencil_ind);
            free(values);

            HYPRE_StructMatrix R_l = pfmg_data->RT_l[space_level];
            stencil = hypre_StructMatrixStencil(R_l);
            grid = hypre_StructMatrixGrid(R_l);
            stencil_size = hypre_StructStencilSize(stencil);
            nrows_l = hypre_StructGridLocalSize(grid);
            nnz_l = nrows_l * stencil_size;
            values = (double *)calloc(nnz_l, sizeof(double));
            ilow = grid->boxes->boxes[0].imin;
            iup = grid->boxes->boxes[0].imax;
            stencil_ind = (int *)calloc(stencil_size, sizeof(int));
            for (int j = 0; j < stencil_size; j++) stencil_ind[j] = j;

            HYPRE_StructMatrixCreate(comm,
                                     grid,
                                     stencil,
                                     &(R_lev[level]));
            HYPRE_StructMatrixInitialize(R_lev[level]);
            HYPRE_StructMatrixGetBoxValues(R_l, ilow, iup, stencil_size, stencil_ind, values);
            for (int i = 0; i < nnz_l; i++) values[i] = .5;
            HYPRE_StructMatrixSetBoxValues(R_lev[level], ilow, iup, stencil_size, stencil_ind, values);
            //hypre_StructMatrixClearBoundary(R_lev[level]);
            HYPRE_StructMatrixAssemble(R_lev[level]);

            free(stencil_ind);
            free(values);
         }

         x_lev[level] = SetupVectorLevel(A_imp_lev[level], level);
         b_lev[level] = SetupVectorLevel(A_imp_lev[level], level);
      }
      
      //for (int level = 0; level < num_levels-1; level++){
      //   char filename[100];
      //   sprintf(filename, "../../../matlab_codes/A_%d", level);
      //   HYPRE_StructMatrixPrint(filename, A_imp_lev[level], 0);
      //   sprintf(filename, "../../../matlab_codes/R_%d", level);
      //   HYPRE_StructMatrixPrint(filename, R_lev[level], 0);
      //   sprintf(filename, "../../../matlab_codes/P_%d", level);
      //   HYPRE_StructMatrixPrint(filename, P_lev[level], 0);
      //}
      //int level = num_levels-1;
      //char filename[100];
      //sprintf(filename, "../../../matlab_codes/A_%d", level);
      //HYPRE_StructMatrixPrint(filename, A_imp_lev[level], 0);
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

      //for (int level = 0; level < num_levels; level++){
      //   char filename[100];
      //   sprintf(filename, "../../../matlab_codes/A_%d", level);
      //   HYPRE_StructMatrixPrint(filename, A_imp_lev[level], 0);
      //}
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
   HYPRE_StructPFMGSetRAPType(*pfmg_solver, 0);
   HYPRE_StructPFMGSetRelaxType(*pfmg_solver, relax_type);
   HYPRE_StructPFMGSetJacobiWeight(*pfmg_solver, jacobi_weight);
   HYPRE_StructPFMGSetNumPreRelax(*pfmg_solver, n_pre);
   HYPRE_StructPFMGSetNumPostRelax(*pfmg_solver, n_post);
   HYPRE_StructPFMGSetSkipRelax(*pfmg_solver, 0);
   HYPRE_StructPFMGSetLogging(*pfmg_solver, 1);
   HYPRE_StructPFMGSetPrintLevel(*pfmg_solver, 1);
}

void HypreSolver::Restrict(HYPRE_StructVector y_f, HYPRE_StructVector y_c, int f_level, int c_level)
{
   //if (hypre_StructMatrixConstantCoefficient(A_imp)){
   //   hypre_StructVectorClearAllValues(y_c);
   //}
   hypre_StructVectorSetConstantValues(y_c, 0.0);
   hypre_SemiRestrict(restrict_data_lev[f_level], R_lev[f_level], y_f, y_c);
   HYPRE_StructVectorScaleValues(y_c, .5);
}

void HypreSolver::Prolong(HYPRE_StructVector y_f, HYPRE_StructVector y_c, int f_level, int c_level)
{
   //if (hypre_StructMatrixConstantCoefficient(A_imp)){
   //   hypre_StructVectorClearAllValues(y_f);
   //}
   hypre_StructVectorSetConstantValues(y_f, 0.0);
   hypre_SemiInterp(interp_data_lev[f_level], P_lev[f_level], y_c, y_f);
   //HYPRE_StructVectorScaleValues(y_f, 2.0); 
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

 
   //if (coarse_space_data != NULL) HYPRE_StructPFMGDestroy(coarse_space_data);
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
