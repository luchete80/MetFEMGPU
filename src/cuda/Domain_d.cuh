#ifndef _DOMAIN_D_CUH_
#define _DOMAIN_D_CUH_

#include <cuda.h>
#include "cudautils.cuh"

namespace MetFEM{
class Domain_d {
public:
	void Domain_d::SetDimension(const int &node_count, const int &elem_count); //ELEM TYPE???
  void Domain_d::AddBoxLength(double3 const & V, double3 const & L, const double &r);
  
  __device__ void calcElemJacobian ();
  __device__ void calcElemJAndDerivatives/*_FullInt*/();
	
	int threadsPerBlock, blocksPerGrid; //TO BE USED BY SOLVER
	
	const int & getElemCount()const{return m_elem_count;}
	const int & getNodeCount()const{return m_node_count;}
	
	void SolveChungHulbert();
	
protected:
	int 						m_dim;
  int             m_node_count, m_elem_count;
	
	unsigned int 		*m_elnod;
	
	double3* 				x; //Vector is double
	double3* 				v;
	double3* 				a;
	double3* 				u;

	double 					*p;
  
  bool            m_red_int;  //Reduced integration, 1 GAUSS POINT
  int             m_gp_count; //Assuming constant gauss points
  int             m_nodxelem;
	
  //Updated lagrangian formulation
  //real(fp_kind), dimension(:,:,:,:), allocatable :: BL,BNL, jacob, dHxy,dHxy_detJ, dHxy0,math, dHrs !!!DIM: e,gp,,:,:, is it necesary to store dHrs??? is only because is used twice, at J and dHxy
  double         *dHxy_detJ ; //
	
	
	double					Time;    				//Current time of simulation at each solving step
	double					deltat;					//Time Step
	double					deltatmin;			//Minimum Time Step
	double					deltatint;			//Initial Time Step
	
	
	
};

__global__ void calcElemJAndDerivKernel(Domain_d *dom_d);

}; //Namespace

#endif