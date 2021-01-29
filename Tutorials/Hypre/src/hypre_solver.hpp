#ifndef HYPRE_SOLVER_HPP
#define HYPRE_SOLVER_HPP

#include "hypre_struct.hpp"
#include "smg.h"
#include "pfmg.h"

#define SOLVER_PFMG 0
#define SOLVER_PFMG_CG 1
#define SOLVER_JACOBI 2
#define SOLVER_JACOBI_CG 3

#define RELAX_JACOBI 0
#define RELAX_WEIGHTED_JACOBI 1
#define RELAX_RBGS 2
#define RELAX_RBGS_NONSYMMETRIC 3

using namespace std;

class HypreSolver: public HypreStruct {
   protected:
      ;

   public:
      HYPRE_StructMatrix A_exp = NULL;
      HYPRE_StructMatrix A_imp = NULL;
      HYPRE_StructVector b;
      HYPRE_StructVector x;
      HYPRE_StructSolver solver_imp = NULL;
      HYPRE_StructSolver precond_imp = NULL;

      hypre_PFMGData *pfmg_data;
      HYPRE_StructSolver pfmg_level_data;
      HYPRE_StructMatrix *A_exp_lev;
      HYPRE_StructMatrix *A_imp_lev;
      hypre_StructMatrix **RT_lev;
      hypre_StructMatrix **P_lev;
      hypre_StructVector **x_lev;
      hypre_StructVector **b_lev;
      void **restrict_data_lev;
      void **interp_data_lev;
      int *nrows_lev;
      int **ilower_lev;
      int **iupper_lev;

      int num_levels;
      int n_pre, n_post;
      int smg_print_level, pcg_print_level;
      double tol;
      int max_iter;
      double jacobi_weight;
      int relax_type;
      int solver_type;
      int max_levels;
      int setup_done = 0;
      //int setup_spacial_solver_in_FComp_flag = 0;

      HypreSolver(MPI_Comm in_comm = MPI_COMM_WORLD, int dim = 2, int max_iter = 100, int max_levels = 0, int num_grid_points = 32);
      ~HypreSolver(void);

      int GetNumRowsLevel(int level_index);
      int GetNumLevels(void);
      void SetupStructSolver(HYPRE_StructMatrix *A, HYPRE_StructSolver *solver, HYPRE_StructSolver *precond);
      void CleanupStructSolver(HYPRE_StructSolver *solver, HYPRE_StructSolver *precond);
      void SetupStructPFMGSolver(HYPRE_StructSolver *pfmg_solver);
      void SetupSpacialCoarsen(void);
      void SetupMatrix(HYPRE_StructMatrix *A,
                       int num_grid_points,
                       int level_index,
                       int spacial_coarsen_flag,
                       int implicit_flag,
                       double dtq);
      double* UpdateImplicitMatrix(HYPRE_StructMatrix *A, double dtq);
      HYPRE_StructVector SetupVector(void);
      void FEval(double *y, double t, int level_index, double **f);
      void FComp(double **y, double t, double dtq, double *rhs, int level_index, double **f);
      void Restrict(HYPRE_StructVector y_f, HYPRE_StructVector y_c, int f_level, int c_level);
      void Prolong(HYPRE_StructVector y_f, HYPRE_StructVector y_c, int f_level, int c_level);
      void Cleanup(void);
};

#endif
