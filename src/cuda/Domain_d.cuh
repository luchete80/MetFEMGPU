#ifndef _DOMAIN_D_CUH_
#define _DOMAIN_D_CUH_

#include <cuda.h>
#include <stdio.h>
#include "cudautils.cuh"
#include "Material.cuh"

class Matrix;


namespace MetFEM{
class Domain_d {
public:
	void Domain_d::SetDimension(const int &node_count, const int &elem_count); //ELEM TYPE???
  void Domain_d::AddBoxLength(double3 const & V, double3 const & L, const double &r);
  
  __device__ void calcElemJacobian ();
  __device__ void calcElemJAndDerivatives/*_FullInt*/();

	__device__ void calcElemStrains();
  __device__ void calcElemForces();
  __device__ void calcElemPressure(); //FROM STRAIN

  __device__ void assemblyForces();
  
  inline __device__ double & getDerivative(const int &e, const int &gp, const int &i, const int &j); //I AND J ARE: DIMENSION AND NODE
  inline __device__ void     setDerivative(const int &e, const int &gp, const int &i, const int &j, const double &); //I AND J ARE: DIMENSION AND NODE
  inline __device__ void     setDerivative(const int &e, const int &gp, Matrix *); //I AND J ARE: DIMENSION AND NODE
  
	int threadsPerBlock, blocksPerGrid; //TO BE USED BY SOLVER
	
	const int & getElemCount()const{return m_elem_count;}
	const int & getNodeCount()const{return m_node_count;}
  
  inline __device__ double & getSigma  (const int &e, const int &gp, const int &i, const int &j){if (j<m_dim) return m_sigma   [e*m_gp_count*6 + symm_idx[i][j]];}
  inline __device__ double & getStrRate(const int &e, const int &gp, const int &i, const int &j){if (j<m_dim) return m_str_rate[e*m_gp_count*6 + symm_idx[i][j]];}
	
  //__device__ double3 & getVElem(const int &e, const int &n){return v[m_elnod[e*m_nodxelem+n]];}
  inline __device__ double  getVElem(const int &e, const int &n,const int &d){return v[m_elnod[n]+d];}  
  
	void SolveChungHulbert();
	
protected:
  int symm_idx[3][3] = {{0,3,4},{3,1,5},{4,5,2}};
	int 						m_dim;
  int             m_node_count, m_elem_count;
	
	unsigned int 		*m_elnod;
  //unsigned int    *m_eloffset; //FROM FIRST ELEMENT NODE
	
	double3* 				x; //Vector is double
	double* 				v; //CHANGED TO DOUBLE
	double3* 				a;
	double3* 				u;

	double 					*p;

  Material_ **mat; //pointer to material of each particle
  Material_ *materials; //All materials 
  
  bool            m_red_int;  //Reduced integration, 1 GAUSS POINT
  int             m_gp_count; //Assuming constant gauss points
  int             m_nodxelem;
  
  /////////////////////// //LOW ACCESS SPEED; BUT NOT DYNAMIC CREATION ///////////////////////////////////
  Matrix          **m_dHrs;     //LOW ACCESS SPEED; BUT NOT DYNAMIC CREATION
  Matrix          **x2;         //LOW ACCESS SPEED; BUT NOT DYNAMIC CREATION
  Matrix          **dH_dxyz; 
  
  Matrix          *m_jacob;
  
  double          *m_str_rate, *m_rot_rate;
  double          *m_f_elem;  //Necesary?
  double          *m_f;       //NODAL
  double          *m_sigma;
	
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

__global__ void assemblyForcesKernel(Domain_d *dom_d);
__global__ void calcElemForcesKernel (Domain_d *);
__global__ void calcElemPressureKernel (Domain_d *);

inline __device__ double & Domain_d::getDerivative(const int &e, const int &gp, const int &i, const int &j){ //I AND J ARE: DIMENSION AND NODE
      if (i == 0)     return m_dH_detJ_dx[e*(m_nodxelem * m_gp_count) + gp * m_gp_count + i];
      else if (i==1)  return m_dH_detJ_dx[e*(m_nodxelem * m_gp_count) + gp * m_gp_count + i];
      else if (i==2)  return m_dH_detJ_dx[e*(m_nodxelem * m_gp_count) + gp * m_gp_count + i];
      else printf ("ERROR: WROWNG DERIVATIVE DIMENSION.");
}

inline __device__ void Domain_d::setDerivative(const int &e, const int &gp, const int &i, const int &j, const double &v){ //I AND J ARE: DIMENSION AND NODE
      if (i == 0)     m_dH_detJ_dx[e*(m_nodxelem * m_gp_count) + gp * m_gp_count + j] = v;
      else if (i==1)  m_dH_detJ_dx[e*(m_nodxelem * m_gp_count) + gp * m_gp_count + j] = v;
      else if (i==2)  m_dH_detJ_dx[e*(m_nodxelem * m_gp_count) + gp * m_gp_count + j] = v;
      else printf ("ERROR: WRONG DERIVATIVE DIMENSION.");
}

inline __device__ void Domain_d::setDerivative(const int &e, const int &gp, Matrix *m){ //I AND J ARE: DIMENSION AND NODE
      // for (int j = 0;j<3;j++)  m_dH_detJ_dx[e*(m_nodxelem * m_gp_count) + gp * m_gp_count + i] = v;
      // else if (i==1)  m_dH_detJ_dx[e*(m_nodxelem * m_gp_count) + gp * m_gp_count + i] = v;
      // else if (i==2)  m_dH_detJ_dx[e*(m_nodxelem * m_gp_count) + gp * m_gp_count + i] = v;
      //else printf ("ERROR: WRONG DERIVATIVE DIMENSION.");
}

}; //Namespace
#endif