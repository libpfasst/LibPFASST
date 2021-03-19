
#include "hypre_solver.hpp"

int main(int argc, char *argv[])
{
   int num_procs, rank;
   MPI_Init(NULL, NULL);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   HypreSolver *hypre_solver;
   int level_index = 0;
   int n = 4;
   int num_grid_points;

   int arg_index = 0;

   while (arg_index < argc){
      if (strcmp(argv[arg_index], "-n") == 0){
         arg_index++;
	 n = atoi(argv[arg_index]);
      }
      arg_index++;
   }

   num_grid_points = n / (int)sqrt((double)num_procs); 

   hypre_solver = new HypreSolver(MPI_COMM_WORLD, 2, 50, 15);
   hypre_solver->SetupMatrix(num_grid_points, level_index, 0); 

   int num_rows = hypre_solver->GetNumRows();
   double t = 1.0, dtq = .1;
   double *y = (double *)malloc(num_rows * sizeof(double));
   double *f = (double *)malloc(num_rows * sizeof(double));
   double *rhs = (double *)malloc(num_rows * sizeof(double));
   for (int i = 0; i < num_rows; i++){
      rhs[i] = rank*num_rows + i;
   }
   //hypre_solver->FComp(&y, t, dtq, rhs, level_index, &f);

   hypre_solver->FEval(rhs, t, level_index, &f);
   MPI_Finalize();
   return 0;

   hypre_solver->Cleanup();

   //for (int i = 0; i < num_rows; i++) printf("%f\n", f[i]);

   delete hypre_solver;

   free(y);
   free(f);
   free(rhs);
   return 0;
}
