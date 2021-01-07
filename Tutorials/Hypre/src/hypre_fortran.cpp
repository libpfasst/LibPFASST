#include "hypre_fortran.hpp"

/*
 *  C++ functions called from Fortran to manipulate Hypre vector and solver objects
 */

HypreSolver *glob_hypre_solver;
int glob_spacial_coarsen_flag;

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
      //if (glob_spacial_coarsen_flag == 1){
         level_index = PfasstToHypreLevelIndex(pfasst_level_index, hypre_solver->GetNumLevels());
      //}
      //else {
      //   level_index = 0;
      //}
      return hypre_solver->GetNumRowsLevel(level_index);
   }

   int HypreSolverGetExtentLevel(HypreSolver *hypre_solver, int pfasst_level_index, int i)
   {
      int level_index;
      //if (glob_spacial_coarsen_flag == 1){
         level_index = PfasstToHypreLevelIndex(pfasst_level_index, hypre_solver->GetNumLevels());
      //}
      //else {
      //   level_index = 0;
      //}
      int extent;
      switch(i){
         case 0:
            extent = hypre_solver->ilower_l[level_index][0];
            break;
         case 1:
            extent = hypre_solver->ilower_l[level_index][1];
            break;
         case 2:
            extent = hypre_solver->iupper_l[level_index][0];
            break;
         case 3:
            extent = hypre_solver->iupper_l[level_index][1];
            break;
      }
      return extent;
   }

   void HypreVectorCreate(HypreVector **hypre_vector, 
                          int num_grid_points,
                          int comm_color,
                          int space_dim,
                          int nrows,
                          int ilower0,
                          int ilower1,
                          int iupper0,
                          int iupper1)
   {
      int extents[4] = {ilower0, ilower1, iupper0, iupper1};
      MPI_Comm newcomm;
      MPI_Comm_split(MPI_COMM_WORLD, comm_color, 0, &newcomm);
      *hypre_vector = new HypreVector(num_grid_points, 0.0, newcomm, space_dim, nrows, extents);
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

   void HypreSolverInit(HypreSolver **hypre_solver, int pfasst_level_index, int num_grid_points, int comm_color, int space_dim, int max_iter, int max_levels, int spacial_coarsen_flag)
   {
      int level_index = PfasstToHypreLevelIndex(pfasst_level_index, max_levels);
      MPI_Comm newcomm;
      MPI_Comm_split(MPI_COMM_WORLD, comm_color, 0, &newcomm);
      *hypre_solver = new HypreSolver(newcomm, space_dim, max_iter, max_levels);
      (*hypre_solver)->SetupMatrix(num_grid_points, level_index, spacial_coarsen_flag);
      if (level_index == 0){
         glob_hypre_solver = *hypre_solver;
         glob_spacial_coarsen_flag = spacial_coarsen_flag;
      }
   }

   void HypreSolverDestroy(HypreSolver *hypre_solver, int pfasst_level_index)
   {
      int level_index = PfasstToHypreLevelIndex(pfasst_level_index, hypre_solver->GetNumLevels());
      delete hypre_solver;
   }

   void HypreSolverFEval(HypreSolver *hypre_solver, HypreVector *y, double t, int pfasst_level_index, HypreVector *f, int piece)
   {
      int level_index = PfasstToHypreLevelIndex(pfasst_level_index, hypre_solver->GetNumLevels());

      int nrows = hypre_solver->GetNumRows();
      double *f_values = (double *)malloc(nrows * sizeof(double));
      if (piece == 1){
         for (int i = 0; i < nrows; i++){
            f_values[i] = 0.0;
         }
      }
      else if (piece == 2){
         hypre_solver->FEval(y->GetBoxValues(), t, level_index, &f_values);
      }
      else {
         cout << "ERROR in \"HypreSolverFEval()\": Bad value for variable 'piece'\n";
         exit(1);
      }
      f->SetBoxValues(f_values);
   }

   void HypreSolverFComp(HypreSolver *hypre_solver, HypreVector *y, double t, double dtq, HypreVector *rhs, int pfasst_level_index, HypreVector *f, int piece)
   {
      if (piece != 2){
         cout << "ERROR in \"HypreSolverFComp()\": Bad value for variable 'piece'\n";
         exit(1);
      }

      int level_index = PfasstToHypreLevelIndex(pfasst_level_index, hypre_solver->GetNumLevels());

      int nrows = hypre_solver->GetNumRows();
      double *y_values = (double *)malloc(nrows * sizeof(double));
      double *f_values = (double *)malloc(nrows * sizeof(double));
      hypre_solver->FComp(&y_values, t, dtq, rhs->GetBoxValues(), level_index, &f_values);
      y->SetBoxValues(y_values);
      f->SetBoxValues(f_values);
   }

   double HypreMaxErr(HypreVector *hypre_vector, double t, double init_cond)
   {
      int Q = 100, P = 100;
      
      double *u = hypre_vector->HeatEquTrueSol(t, P, Q, init_cond);
      double *x = hypre_vector->GetBoxValues();
      
      double max_err = 0.0;
      for (int i = 0; i < hypre_vector->GetNumRows(); i++){
         double abs_diff = abs(x[i] - u[i]);
         //printf("%f %f\n", x[i], u[i]);
         if (abs_diff > max_err){
            max_err = abs_diff;
         }
      }

      free(u);
 
      return max_err;
   }

   void HypreRestrict(HypreVector *y_f, HypreVector *y_c, int pfasst_f_level, int pfasst_c_level)
   {
      if (glob_spacial_coarsen_flag == 1){   
         int f_level = PfasstToHypreLevelIndex(pfasst_f_level, glob_hypre_solver->GetNumLevels());
         int c_level = PfasstToHypreLevelIndex(pfasst_c_level, glob_hypre_solver->GetNumLevels());
         glob_hypre_solver->Restrict(y_f->GetVector(), y_c->GetVector(), f_level, c_level);
      }
      else {
         HypreVectorCopy(y_c, y_f);
      }
   }

   void HypreProlong(HypreVector *y_f, HypreVector *y_c, int pfasst_f_level, int pfasst_c_level)
   {
      if (glob_spacial_coarsen_flag == 1){
         int f_level = PfasstToHypreLevelIndex(pfasst_f_level, glob_hypre_solver->GetNumLevels());
         int c_level = PfasstToHypreLevelIndex(pfasst_c_level, glob_hypre_solver->GetNumLevels());
         glob_hypre_solver->Prolong(y_f->GetVector(), y_c->GetVector(), f_level, c_level);
      }
      else {
         HypreVectorCopy(y_f, y_c);
      }
   }

   void HypreSolverSetLevelData(HypreSolver **y, HypreSolver *x, int pfasst_level_index)
   {
      int level_index;
      if (glob_spacial_coarsen_flag == 1){
         level_index = PfasstToHypreLevelIndex(pfasst_level_index, x->GetNumLevels());
         (*y)->A = x->A_l[level_index];
         (*y)->x = x->x_l[level_index];
         (*y)->b = x->b_l[level_index];
         (*y)->nrows = x->nrows_l[level_index];
         (*y)->ilower = x->ilower_l[level_index];
         (*y)->iupper = x->iupper_l[level_index];
         (*y)->nentries = x->nentries;
         (*y)->stencil_indices = x->stencil_indices;
         (*y)->stencil = hypre_StructMatrixStencil((*y)->A);
         (*y)->nnz = (*y)->nrows * (*y)->nentries;
         (*y)->num_levels = x->num_levels;
      }
   }

   int HypreSolverGetNumLevels(HypreSolver *hypre_solver)
   {
      return hypre_solver->GetNumLevels();
   }
}
