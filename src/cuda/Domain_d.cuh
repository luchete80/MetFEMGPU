#ifndef _DOMAIN_D_CUH_
#define _DOMAIN_CUH_

#include <cuda.h>

namespace MetFEM{
class Domain_d {
public:
	void Domain_d::SetDimension(const int &node_count, const int &elem_count); //ELEM TYPE???

protected:
	int 						m_dim;
	
	double3* 				x; //Vector is double
	double3* 				v;
	double3* 				a;
	double3* 				u;

	double 					*p;
	
	
	
	double					Time;    				//Current time of simulation at each solving step
	double					deltat;					//Time Step
	double					deltatmin;			//Minimum Time Step
	double					deltatint;			//Initial Time Step
	
};

}; //Namespace

#endif