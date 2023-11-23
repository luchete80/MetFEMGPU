#ifndef _DOMAIN_D_CUH_
#define _DOMAIN_CUH_

#include <cuda.h>
#include "cudautils.cuh"

namespace MetFEM{
class Domain_d {
public:
	void Domain_d::SetDimension(const int &node_count, const int &elem_count); //ELEM TYPE???
  void Domain_d::AddBoxLength(double3 const & V, double3 const & L, const double &r);

protected:
	int 						m_dim;
  int             m_node_count, m_elem_count;
	
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