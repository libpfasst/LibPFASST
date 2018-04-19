#!/bin/bash
for value in {2..9}
	 do
	     echo "IMEX SDC with nnodes=" $value
	     ./main.exe sdc.nml nnodes=$value | grep 'step: 008'
	     echo "Implicit SDC with nnodes=" $value
	     ./main.exe sdc.nml nnodes=$value imex_stat=1 | grep 'step: 008'
	     echo "Explicit SDC with nnodes=" $value
	     ./main.exe sdc.nml nnodes=$value imex_stat=0 nsteps=64 nu=0.001 | grep 'step: 032'
done

echo "Imex MLSDC with nnodes=[3 5 9]" 
./main.exe mlsdc.nml  imex_stat=1  | grep 'step: 008'
echo "Implicit MLSDC with nnodes=[3 5 9]" $value
./main.exe mlsdc.nml  imex_stat=2 nnodes=[2 3 5] | grep 'step: 008'



