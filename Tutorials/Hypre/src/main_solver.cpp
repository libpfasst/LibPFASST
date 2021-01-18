#include "hypre_solver.hpp"

int main(void)
{
   HypreSolver *hypre_solver;
   int level_index = 0;
   int num_grid_points = 3;
   
   MPI_Init(NULL, NULL);

   hypre_solver = new HypreSolver;

   hypre_solver->SetComm(MPI_COMM_WORLD);
   hypre_solver->SetupMatrix(level_index, num_grid_points); 

   int n = hypre_solver->GetNumRows();
   int piece = 2;
   double t = 1.0, dtq = .1;
   double *y = (double *)malloc(n * sizeof(double));
   double *f = (double *)malloc(n * sizeof(double));
   double *rhs = (double *)malloc(n * sizeof(double));
   for (int i = 0; i < n; i++){
      rhs[i] = 1.0;
   }
   hypre_solver->FComp(&y, t, dtq, rhs, level_index, &f, piece);

   //hypre_solver->FEval(rhs, t, level_index, &f, piece);

   hypre_solver->Cleanup();

   for (int i = 0; i < n; i++) printf("%f\n", f[i]);

   delete hypre_solver;

   free(y);
   free(f);
   free(rhs);
   MPI_Finalize();
   return 0;
}
