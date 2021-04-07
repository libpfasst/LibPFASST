#ifndef HYPRE_FORTRAN_HPP
#define HYPRE_FORTRAN_HPP

#include "hypre_vector.hpp"
#include "hypre_solver.hpp"
#include <new>

/*
 *  Prototypes for C++ functions that are called from Fortran
 */

extern "C"
{
   int PfasstToHypreLevelIndex(int pfasst_level_index, int num_levels);
   int HypreSolverGetNumRowsLevel(HypreSolver *hypre_solver, int pfasst_level_index);
   void HypreVectorCreate(HypreVector **hypre_vector,
                          int num_grid_points,
                          int comm_color,
                          int space_dim,
                          int nrows,
                          int ilower0,
                          int ilower1,
                          int iupper0,
                          int iupper1);
   void HypreVectorDestroy(HypreVector *hypre_vector);
   void HypreVectorSetVal(HypreVector *hypre_vector, double val);
   void HypreVectorSetSinInitCond(HypreVector *hypre_vector);
   void HypreVectorCopy(HypreVector *dest, HypreVector *src);
   double *HypreVectorPack(HypreVector *hypre_vector);
   void HypreVectorUnpack(HypreVector *hypre_vector, double *z);
   double HypreVectorNorm(HypreVector *hypre_vector);
   void HypreVectorAxpy(HypreVector *y, double a, HypreVector *x);
   void HypreVectorPrint(HypreVector *hypre_vector);
   void HypreSolverInit(HypreSolver **hypre_solver,
                        int pfasst_level_index,
                        int num_grid_points,
                        int comm_color,
                        int space_dim,
                        int max_iter,
                        int max_levels,
                        int spacial_coarsen_flag);
   void HypreImplicitSolverInit(HypreSolver **hypre_solver,
                                int pfasst_level_index,
                                int num_grid_points,
                                int comm_color,
                                int space_dim,
                                int max_iter,
                                int max_levels,
                                double dtq);
   void HypreSolverDestroy(HypreSolver *hypre_solver, int pfasst_level_index);
   void HypreSolverFEval(HypreSolver *hypre_solver, HypreVector *y, double t, int pfasst_level_index, HypreVector *f, int piece);
   void HypreSolverFComp(HypreSolver *hypre_solver, HypreVector *y, double t, double dtq, HypreVector *rhs, int pfasst_level_index, HypreVector *f, int piece);
   double HypreMaxErr(HypreVector *hypre_vector, double t, double init_cond);
   void HypreRestrict(HypreVector *y_f, HypreVector *y_c, int f_level, int c_level);
   void HypreProlong(HypreVector *y_f, HypreVector *y_c, int f_level, int c_level);
   int HypreSolverGetExtentLevel(HypreSolver *hypre_solver, int pfasst_level_index, int i);
   void HypreSolverSetLevelData(HypreSolver **y, HypreSolver *x, int pfasst_level_index);
   int HypreSolverGetNumLevels(HypreSolver *hypre_solver);
   void GetHypreStats(void);
}

#endif
