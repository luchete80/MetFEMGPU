#ifndef _DOMAIN_D_CUH_
#define _DOMAIN_D_CUH_

#include <cuda.h>
#include "cudautils.cuh"

class Matrix;

namespace MetFEM{
class Domain_d {
public:
	void Domain_d::SetDimension(const int &node_count, const int &elem_count); //ELEM TYPE???
  void Domain_d::AddBoxLength(double3 const & V, double3 const & L, const double &r);
  
  __device__ void calcElemJacobian ();
  __device__ void calcElemJAndDerivatives/*_FullInt*/();

	__device__ void calcElemStrains();
  
  __device__ double & getDerivative(const int &e, const int &gp, const int &i, const int &j); //I AND J ARE: DIMENSION AND NODE
  
	int threadsPerBlock, blocksPerGrid; //TO BE USED BY SOLVER
	
	const int & getElemCount()const{return m_elem_count;}
	const int & getNodeCount()const{return m_node_count;}
	
  __device__ double & getVElem(const int &e, const int &n){return v[]}
  
  
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
  
  /////////////////////// //LOW ACCESS SPEED; BUT NOT DYNAMIC CREATION ///////////////////////////////////
  Matrix          **m_dHrs;     //LOW ACCESS SPEED; BUT NOT DYNAMIC CREATION
  Matrix          **x2;         //LOW ACCESS SPEED; BUT NOT DYNAMIC CREATION
  Matrix          **dH_dxyz; 
  
  Matrix          *m_jacob;
  
  double          *m_strain_rate;
	
  //Updated lagrangian formulation
  //real(fp_kind), dimension(:,:,:,:), allocatable :: BL,BNL, jacob, dHxy,dHxy_detJ, dHxy0,math, dHrs !!!DIM: e,gp,,:,:, is it necesary to store dHrs??? is only because is used twice, at J and dHxy
  double         *dHxy_detJ ; //
  double         *m_detJ;
  
  ////// THESE ARE SUBDIVIDED FOR HIGHER ACCESS SPEED (LOWER OFFSET)
  double         *m_dH_detJ_dx ; //GAUSS_POINT, AGRUPATED BY ELEM 1 GP1, ELEM 1 GP 2 .... ELEM 2 GP 1..
  double         *m_dH_detJ_dy ; 	
  double         *m_dH_detJ_dz ; 
	
	double					Time;    				//Current time of simulation at each solving step
	double					deltat;					//Time Step
	double					deltatmin;			//Minimum Time Step
	double					deltatint;			//Initial Time Step
	
	
	
};

__global__ void calcElemJAndDerivKernel(Domain_d *dom_d);
__global__ void calcElemStrainsKernel(Domain_d *dom_d);

}; //Namespace

#endif