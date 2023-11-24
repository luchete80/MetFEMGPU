#ifndef _DOMAIN_D_CUH_
#define _DOMAIN_CUH_

#include <cuda.h>
#include "cudautils.cuh"

namespace MetFEM{
class Domain_d {
public:
	void Domain_d::SetDimension(const int &node_count, const int &elem_count); //ELEM TYPE???
  void Domain_d::AddBoxLength(double3 const & V, double3 const & L, const double &r);
  
  __device__ void calcElemJacobian ();
  __device__ void calcDerivatives_FullInt();

protected:
	int 						m_dim;
  int             m_node_count, m_elem_count;
	
	unsigned int 		*elnod;
	
	double3* 				x; //Vector is double
	double3* 				v;
	double3* 				a;
	double3* 				u;

	double 					*p;
  
  bool            red_int;  //Reduced integration, 1 GAUSS POINT
  int             gp_count; //Assuming constant gauss points
  int             nodxelem;
	
  //Updated lagrangian formulation
  //real(fp_kind), dimension(:,:,:,:), allocatable :: BL,BNL, jacob, dHxy,dHxy_detJ, dHxy0,math, dHrs !!!DIM: e,gp,,:,:, is it necesary to store dHrs??? is only because is used twice, at J and dHxy
  double         *dHxy_detJ ; //
	
	
	double					Time;    				//Current time of simulation at each solving step
	double					deltat;					//Time Step
	double					deltatmin;			//Minimum Time Step
	double					deltatint;			//Initial Time Step
	
};

}; //Namespace

#endif