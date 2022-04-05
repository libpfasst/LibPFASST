#ifndef HYPRE_VECTOR_HPP
#define HYPRE_VECTOR_HPP

#include "hypre_struct.hpp"

using namespace std;

class HypreVector: public HypreStruct {
   protected:
      HYPRE_StructVector v;
   public:
      HypreVector(int num_grid_points = 10, double value = 0.0, MPI_Comm in_comm = MPI_COMM_WORLD, int in_dim = 2, int in_nrows = -1, int *extents = NULL);
      ~HypreVector(void);

      HYPRE_StructVector GetVector(void);
      double *GetBoxValues(void);
      void SetBoxValues(double *values);
      void SetVal(double val);
      void SetSinInitCond(void);
      double Norm(void);
      void Axpy(double alpha, HypreVector *x);
      void Print(void);
};

#endif
