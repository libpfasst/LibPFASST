#include "hypre_fortran.hpp"
#include "hypre_solver.hpp"

/*
 *  C++ functions called from Fortran to manipulate Hypre vector and solver objects
 */

HypreSolver *glob_hypre_solver = nullptr;
int glob_spatial_coarsen_flag = 0;
MPI_Comm glob_space_comm = MPI_COMM_NULL;

int FComp_count = 0;
double FComp_wtime = 0;

/* TODO: handle case when number of pfasst levels > number of hypre levels */
extern "C"
{
   int PfasstToHypreLevelIndex(int pfasst_level_index, int num_levels)
   {
      return abs(num_levels - pfasst_level_index);
   }

   int HypreSolverGetNumRowsLevel(HypreSolver *hypre_solver, int pfasst_level_index)
   {
      int level_index;
      int num_rows;

      level_index = PfasstToHypreLevelIndex(pfasst_level_index, glob_hypre_solver->GetNumLevels());
      num_rows = glob_hypre_solver->GetNumRowsLevel(level_index);

      return num_rows;
   }

   int HypreSolverGetExtentLevel(HypreSolver *hypre_solver, int pfasst_level_index, int i)
   {
      int level_index, extent;
      level_index = PfasstToHypreLevelIndex(pfasst_level_index, glob_hypre_solver->GetNumLevels());
      switch(i){
         case 0:
            extent = glob_hypre_solver->ilower_lev[level_index][0];
            break;
         case 1:
            extent = glob_hypre_solver->ilower_lev[level_index][1];
            break;
         case 2:
            extent = glob_hypre_solver->iupper_lev[level_index][0];
            break;
         case 3:
            extent = glob_hypre_solver->iupper_lev[level_index][1];
            break;
      }

      return extent;
   }

   void HypreVectorCreate(HypreVector **hypre_vector, 
                          int pfasst_level_index,        
                          int nx,
                          int comm_color,
                          int space_dim)
   {
      if (glob_hypre_solver == nullptr){
         printf("ERROR: must set up hypre matrices/solver before vectors.");
         MPI_Finalize();
         exit(1);
      }

      if (*hypre_vector == nullptr){
         int level_index;
         level_index = PfasstToHypreLevelIndex(pfasst_level_index, glob_hypre_solver->GetNumLevels());
         int extents[4] = {glob_hypre_solver->ilower_lev[level_index][0], 
                           glob_hypre_solver->ilower_lev[level_index][1],
                           glob_hypre_solver->iupper_lev[level_index][0],
                           glob_hypre_solver->iupper_lev[level_index][1]};
         *hypre_vector = new HypreVector(nx,
                                         0.0,
                                         glob_space_comm,
                                         glob_hypre_solver->dim,
                                         glob_hypre_solver->nrows_lev[level_index],
                                         extents,
                                         hypre_StructMatrixGrid(glob_hypre_solver->A_imp_lev[level_index]),
                                         hypre_StructMatrixStencil(glob_hypre_solver->A_imp_lev[level_index]),
                                         0);
      }
   }

   void HypreVectorDestroy(HypreVector *hypre_vector)
   {
      if (hypre_vector != nullptr){
         delete hypre_vector;
      }
   }

   void HypreVectorSetVal(HypreVector *hypre_vector, double val)
   {
      hypre_vector->SetVal(val);
   }

   void HypreVectorSetSinInitCond(HypreVector *hypre_vector)
   {
      hypre_vector->SetSinInitCond();
   }

   void HypreVectorCopy(HypreVector *dest, HypreVector *src)
   {
      HYPRE_StructVectorCopy(src->GetVector(), dest->GetVector());
   }

   double *HypreVectorPack(HypreVector *hypre_vector)
   {
      return hypre_vector->GetBoxValues();
   }

   void HypreVectorUnpack(HypreVector *hypre_vector, double *z)
   {
      hypre_vector->SetBoxValues(z);
   }

   double HypreVectorNorm(HypreVector *hypre_vector)
   {
      return hypre_vector->Norm();
   }

   void HypreVectorAxpy(HypreVector *y, double a, HypreVector *x)
   {
      y->Axpy(a, x);
   }

   void HypreVectorPrint(HypreVector *hypre_vector)
   {
      hypre_vector->Print();
   }

   void HypreSolverInit(HypreSolver **hypre_solver,
                        int pfasst_level_index,
                        int nx,
                        int comm_color,
                        int space_dim,
                        int max_iter,
                        int max_levels,
                        int spatial_coarsen_flag)
   {
      int level_index = PfasstToHypreLevelIndex(pfasst_level_index, max_levels);
      if (glob_space_comm == MPI_COMM_NULL){
         MPI_Comm_split(MPI_COMM_WORLD, comm_color, 0, &glob_space_comm);
      }
      MPI_Comm newcomm = glob_space_comm;

      double cx = 0.0, cy = 0.0;
      int adv_flag = 0;

      if (glob_hypre_solver == nullptr){
         glob_spatial_coarsen_flag = spatial_coarsen_flag;
         glob_hypre_solver = new HypreSolver(newcomm, space_dim, max_iter, max_levels, nx, cx, cy, adv_flag);
         glob_hypre_solver->SetupMatrix(&(glob_hypre_solver->A_exp), nx, level_index, spatial_coarsen_flag, 0, 0.0);
         glob_hypre_solver->SetupMatrix(&(glob_hypre_solver->A_imp), nx, level_index, spatial_coarsen_flag, 0, 0.0);
         glob_hypre_solver->x = glob_hypre_solver->SetupVector();
         glob_hypre_solver->b = glob_hypre_solver->SetupVector();
         glob_hypre_solver->SetupLevels(spatial_coarsen_flag, nx);
      }
      
      if (*hypre_solver == nullptr){
         *hypre_solver = new HypreSolver(newcomm, space_dim, max_iter, max_levels, nx, 0.0, 0.0, adv_flag);
      }
   }

   void HypreImplicitSolverInit(HypreSolver **hypre_solver,
                                int pfasst_level_index,
                                int nx,
                                int comm_color,
                                int space_dim,
                                int max_iter,
                                int max_levels,
                                double dtq)
   {
      int level_index = PfasstToHypreLevelIndex(pfasst_level_index, max_levels);
      //if (glob_spatial_coarsen_flag == 1){
      //   if (level_index >= glob_hypre_solver->spatial_num_levels-1) return;
      //}

      double *dummy = glob_hypre_solver->UpdateImplicitMatrix(level_index, dtq);
      free(dummy);
      glob_hypre_solver->SetupStructSolver(level_index);
   }

   void HypreSolverDestroy(HypreSolver *hypre_solver, int pfasst_level_index)
   {
      int level_index = PfasstToHypreLevelIndex(pfasst_level_index, glob_hypre_solver->GetNumLevels());
      if (hypre_solver != nullptr){
         delete hypre_solver;
      }
   }

   void HypreSolverFEval(HypreSolver *hypre_solver, HypreVector *y, double t, int pfasst_level_index, HypreVector *f, int piece)
   {
      int level_index = PfasstToHypreLevelIndex(pfasst_level_index, glob_hypre_solver->GetNumLevels());
      int nrows = glob_hypre_solver->GetNumRowsLevel(level_index);
      double *f_values = (double *)malloc(nrows * sizeof(double));
      glob_hypre_solver->FEval(y->GetBoxValues(), t, level_index, &f_values, piece);
      f->SetBoxValues(f_values);
      free(f_values);
   }

   void HypreSolverFComp(HypreSolver *hypre_solver, HypreVector *y, double t, double dtq, HypreVector *rhs, int pfasst_level_index, HypreVector *f, int piece)
   {
      if (piece != 2){
         cout << "ERROR in \"HypreSolverFComp()\": Bad value for variable 'piece'\n";
         exit(1);
      }
      
      double wtime_start = MPI_Wtime();

      int level_index = PfasstToHypreLevelIndex(pfasst_level_index, glob_hypre_solver->GetNumLevels());

      int nrows = glob_hypre_solver->GetNumRowsLevel(level_index);
      double *y_values = (double *)malloc(nrows * sizeof(double));
      double *f_values = (double *)malloc(nrows * sizeof(double));
      glob_hypre_solver->FComp(&y_values, t, dtq, rhs->GetBoxValues(), level_index, &f_values);
      y->SetBoxValues(y_values);
      f->SetBoxValues(f_values);
      free(y_values);
      free(f_values);

      //FComp_wtime += MPI_Wtime() - wtime_start;
      FComp_count++;
   }

   double HypreMaxErr(HypreVector *hypre_vector, double t, double init_cond)
   {
      int Q = 100, P = 100;

      double *u = hypre_vector->HeatEquTrueSol(t, P, Q, init_cond);
      double *x = hypre_vector->GetBoxValues();

      double error;
      double max_error = 0.0;
      double error_inner_prod = 0.0, error_L2norm;
      for (int i = 0; i < hypre_vector->GetNumRows(); i++){
         double abs_diff = abs(x[i] - u[i]);
         if (abs_diff > max_error){
            max_error = abs_diff;
         }
         //error_inner_prod += pow(abs_diff, 2.0);
      }

      free(u);

      error = max_error;

      //MPI_Allreduce(&error_inner_prod, &error_L2norm, 1, MPI_DOUBLE, MPI_SUM, hypre_vector->comm);
      //error = sqrt(error_L2norm);

      return error;
   }

   void HypreRestrict(HypreVector *y_f, HypreVector *y_c, int pfasst_f_level, int pfasst_c_level)
   {
      if (glob_spatial_coarsen_flag == 1){   
         int c_level = PfasstToHypreLevelIndex(pfasst_c_level, glob_hypre_solver->GetNumLevels());
         //if (c_level < glob_hypre_solver->spatial_num_levels){
            int f_level = PfasstToHypreLevelIndex(pfasst_f_level, glob_hypre_solver->GetNumLevels());
            glob_hypre_solver->Restrict(y_f->GetVector(), y_c->GetVector(), f_level, c_level);
         //}
         //else {
         //   HypreVectorCopy(y_c, y_f);
         //}
      }
      else {
         HypreVectorCopy(y_c, y_f);
      }
   }

   void HypreProlong(HypreVector *y_f, HypreVector *y_c, int pfasst_f_level, int pfasst_c_level)
   {
      if (glob_spatial_coarsen_flag == 1){
         int c_level = PfasstToHypreLevelIndex(pfasst_c_level, glob_hypre_solver->GetNumLevels());
         //if (c_level < glob_hypre_solver->spatial_num_levels){
            int f_level = PfasstToHypreLevelIndex(pfasst_f_level, glob_hypre_solver->GetNumLevels());
            glob_hypre_solver->Prolong(y_f->GetVector(), y_c->GetVector(), f_level, c_level);
         //}
         //else {
         //   HypreVectorCopy(y_f, y_c);
         //}
      }
      else {
         HypreVectorCopy(y_f, y_c);
      }
   }

   int HypreSolverGetNumLevels(HypreSolver *hypre_solver)
   {
      return glob_hypre_solver->GetNumLevels();
   }

   void GetHypreStats(void)
   {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      printf("rank %d: FComp count %d, wtime %e\n", rank, FComp_count, FComp_wtime);
   }
}
