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

using namespace std;

class HypreSolver: public HypreStruct {
   protected:
      ;

   public:
      /* TODO: make protected */
      HYPRE_StructMatrix A = NULL;
      HYPRE_StructVector b;
      HYPRE_StructVector x;
      HYPRE_StructSolver pfmg_level_data;
      HYPRE_StructSolver solver;
      HYPRE_StructSolver precond;

      hypre_PFMGData *pfmg_data;
      hypre_StructMatrix **A_l;
      hypre_StructMatrix **RT_l;
      hypre_StructMatrix **P_l;
      hypre_StructVector **x_l;
      hypre_StructVector **b_l;
      void **restrict_data_l;
      void **interp_data_l;
      int *nrows_l;
      int **ilower_l;
      int **iupper_l;

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
      /* end TODO */


      HypreSolver(MPI_Comm in_comm = MPI_COMM_WORLD, int dim = 2, int max_iter = 100, int max_levels = 0);
      ~HypreSolver(void);

      HYPRE_StructSolver GetStructSolver(void);
      int GetNumRowsLevel(int level_index);
      int GetNumLevels(void);
      void SetupStructPFMGSolver(HYPRE_StructSolver *pfmg_solver);
      void SetupMatrix(int num_grid_points, int level_index);
      void FEval(double *y, double t, int level_index, double **f);
      void FComp(double **y, double t, double dtq, double *rhs, int level_index, double **f);
      void Restrict(HYPRE_StructVector y_f, HYPRE_StructVector y_c, int f_level, int c_level);
      void Prolong(HYPRE_StructVector y_f, HYPRE_StructVector y_c, int f_level, int c_level);
      void Cleanup(void);
};

#endif
