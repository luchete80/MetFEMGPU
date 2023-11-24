#ifndef _SOLVER_CUH_
#define _SOLVER_CUH_

#include <cuda.h>
#include "cudautils.cuh"

class Domain_d;

namespace MetFEM{

class Solver {
public:
	Solver(){}
	
	virtual Solve(){}
	virtual Solve(Domain_d *dom){}
	
protected:
	Domain_d 	*dom;
};

}; //Namespace

#endif