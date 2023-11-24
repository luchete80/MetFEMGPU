#include "Domain_d.cuh"
#include <iostream>
#include <vector>

#include "tensor.cuh"
#include "Matrix.cuh"

using namespace std;

namespace MetFEM {

void Domain_d::SetDimension(const int &node_count, const int &elem_count){
  
  m_node_count = node_count;
  m_elem_count = elem_count;
  
	cudaMalloc((void **)&x, node_count * sizeof (double3));
	cudaMalloc((void **)&v, node_count * sizeof (double3));
	cudaMalloc((void **)&a, node_count * sizeof (double3));
	
	report_gpu_mem_();

}

void Domain_d::AddBoxLength(double3 const & V, double3 const & L, const double &r){
	    // integer, intent(in):: tag
    // logical, intent(in) :: redint
    // !real(fp_kind), intent(in), allocatable :: V
    // real(fp_kind), dimension(1:3), intent(in)  :: V ! input
    // real(fp_kind), intent(in):: r, Lx, Ly, Lz, Density, h  
    double3 Xp;
    int p, nnodz;
    int nodxelem;
    int nel[3];
    m_dim = 3;
    if (L.z > 0.0) m_dim = 2;
    
    
    nel[0] = (int)(L.x/(2.0*r));
    nel[1] = (int)(L.y/(2.0*r));
    cout << "Nel x: "<<nel[0]<<", y "<<nel[1]<<endl;
    if (m_dim == 2){
      nel[2] = 1;
      nodxelem = 4;
    } else {
      nel[2] = (int)(L.z/(2.0*r));
      nodxelem = 8;
    }
    

    Xp.z = V.z ;
    

    // write (*,*) "Creating Mesh ...", "Elements ", neL.y, ", ",neL.z
  int nc = (nel[0] +1) * (nel[1]+1) * (nel[2]+1);
  int ne = nel[0]*nel[1]*nel[2];
  //thisAllocateNodes((nel[0] +1) * (nel[1]+1) * (nel[2]+1));
    // print *, "Element count in XYZ: ", nel(:)
    // write (*,*) "Box Node count ", node_count

	this->SetDimension(nc,ne);	 //AFTER CREATING DOMAIN
  cout << "Mesh generated. Node count: " << nc<<". Element count: "<<ne<<endl;
  //SPH::Domain	dom;
	//double3 *x =  (double3 *)malloc(dom.Particles.size());
	double3 *x_H =  new double3 [m_node_count];


	//int size = dom.Particles.size() * sizeof(double3);
	cout << "Copying to device..."<<endl;
    
    cout << "Box Particle Count is " << m_node_count <<endl;
    p = 0;
    for (int j = 0; j < (nel[1] +1);j++) {
      Xp.y = V.y;
      for (int j = 0; j < (nel[1] +1);j++){
        Xp.x = V.x;
        for (int i = 0; i < (nel[0] +1);i++){
					//m_node.push_back(new Node(Xp));
					x_H[p] = Xp;
          //nod%x(p,:) = Xp(:);
          cout << "node " << p <<"X: "<<Xp.x<<"Y: "<<Xp.y<<"Z: "<<Xp.z<<endl;
          p++;
          Xp.x = Xp.x + 2.0 * r;
        }
        Xp.y = Xp.y + 2.0 * r;
      }// 
      Xp.z = Xp.z + 2 * r;

    //cout <<"m_node size"<<m_node.size()<<endl;
    } 
		cudaMemcpy(this->x, x_H, m_node_count, cudaMemcpyHostToDevice);    

    // !! ALLOCATE ELEMENTS
    // !! DIMENSION = 2
    int gp = 1;
    if (m_dim == 2) {
      // if (redint .eqv. .False.) then
        // gp = 4
      // end if 
      //call AllocateElements(neL.y * neL.z,gp) !!!!REDUCED INTEGRATION
    } else {
      // if (redint .eqv. .False.) then
        // gp = 8
      // end if 
      // call AllocateElements(neL.y * neL.z*nel(3),gp) 
    }

		unsigned int *elnod_h = new unsigned int [m_elem_count * nodxelem]; //Flattened
    
		int ex, ey, ez;
		std::vector <int> n;
    if (m_dim == 2) {
			n.resize(4);
      int ei = 0;
      for (int ey = 0; ey < nel[1];ey++){
        for (int ex = 0; ex < nel[0];ex++){
        int iv[4];
        elnod_h[ei] = (nel[0]+1)*ey + ex;        iv[ei+1] = (nel[0]+1)*ey + ex+1;
        iv[2] = (nel[0]+1)*(ey+1) + ex+1;        iv[3] = (nel[0]+1)*(ey+1) + ex;
        // cout << i[]
						// n[0]= m_node[iv[0]];
						// n[1]= m_node[(nel[0]+1)*ey + ex+1];
						// n[2]= m_node[(nel[0]+1)*(ey+1)+ex+1];
						// n[3]= m_node[(nel[0]+1)*(ey+1)+ex];
            cout << "Nel x : "<<nel[0]<<endl;
           cout << "nodes "<<endl;
           for (int i=0;i<4;i++)cout << iv[i]<<", ";
						 //m_element.push_back(new El4N2DPE(n));
																							// m_node[(nel[0]+1)*ey + ex+1],
																							// m_node[(nel[0]+1)*(ey+1)+ex+1],
																							// m_node[(nel[0]+1)*(ey+1)+ex]
																							// );
              //elem%elnod(i,:)=[(neL.y+1)*ey + ex+1,(neL.y+1)*ey + ex+2,(neL.y+1)*(ey+1)+ex+2,(neL.y+1)*(ey+1)+ex+1]         
              //print *, "Element ", i , "Elnod", elem%elnod(i,:) 
					 }
      } 
    } else { //dim: 3
      int ei = 0;
      int nnodz = (nel[0]+1)*(nel[1]+1);
      for (int ez = 0; ez < nel[2];ez++)
      for (int ey = 0; ey < nel[1];ey++){
        for (int ex = 0; ex < nel[0];ex++){
          
          int iv[8];
          int nb1 = nnodz*ez + (nel[0]+1)*ey + ex;
          int nb2 = nnodz*ez + (nel[0]+1)*(ey+1) + ex;
          elnod_h[ei  ] = nb1;
          elnod_h[ei+1] = nb1+1;
          elnod_h[ei+2] = nb2+1;
          elnod_h[ei+3] = nb2;
          elnod_h[ei+4] = nb1 + nnodz*(ez+1);
          elnod_h[ei+5] = nb1 + nnodz*(ez+1) + 1;
          elnod_h[ei+6] = nb2 + nnodz*(ez+1) + 1;
          elnod_h[ei+7] = nb2 + nnodz*(ez+1);

          // elem%elnod(i,:) = [ nnodz*ez + (nel(1)+1)*ey + ex+1,nnodz*ez + (nel(1)+1)*ey + ex+2, &
                              // nnodz*ez + (nel(1)+1)*(ey+1)+ex+2,nnodz*ez + (nel(1)+1)*(ey+1)+ex+1, &
                              // nnodz*(ez + 1) + (nel(1)+1)*ey + ex+1,nnodz*(ez + 1) + (nel(1)+1)*ey + ex+2, &
                              // nnodz*(ez + 1) + (nel(1)+1)*(ey+1)+ex+2,nnodz*(ez + 1)+ (nel(1)+1)*(ey+1)+ex+1];
        // cout << i[]
						// n[0]= m_node[iv[0]];
						// n[1]= m_node[(nel[0]+1)*ey + ex+1];
						// n[2]= m_node[(nel[0]+1)*(ey+1)+ex+1];
						// n[3]= m_node[(nel[0]+1)*(ey+1)+ex];
            cout << "Nel x : "<<nel[0]<<endl;
           cout << "nodes "<<endl;
           
           for (int i=0;i<4;i++)cout << iv[i]<<", ";
           ei += nodxelem;
						 //m_element.push_back(new El4N2DPE(n));
																							// m_node[(nel[0]+1)*ey + ex+1],
																							// m_node[(nel[0]+1)*(ey+1)+ex+1],
																							// m_node[(nel[0]+1)*(ey+1)+ex]
																							// );
              //elem%elnod(i,:)=[(neL.y+1)*ey + ex+1,(neL.y+1)*ey + ex+2,(neL.y+1)*(ey+1)+ex+2,(neL.y+1)*(ey+1)+ex+1]         
              //print *, "Element ", i , "Elnod", elem%elnod(i,:) 
					 }
      } 

		}//if dim 
    
    // call AllocateDomain()
    // i = 1
    // do while ( i <= node_count)
      // nod%is_bcv(i,:) = .false.
      // i = i + 1
    // end do
  
    // ! nod%m(:)   = Density * Lx * Ly * Lz / node_count
    // ! nod%rho(:)   = Density
    // elem%rho_0(:,:) = Density
    // !print *, "Particle mass ", nod%m(2)
    
    // !nod%id(:) = tag
    
    // fext_glob (:,:) = 0.0d0
    
    // elem%e_length(:) = Lx !TODO: CHANGE!
    
    // tot_mass = Density * Lx * Ly * Lz
    // if (dim == 2) then !!!assuming plain strain
      // tot_mass = tot_mass / Lz
    // end if
    // print *, "Total Mass: ", tot_mass
    
    // call SearchNodelem
		
		delete [] elnod_h;
}

__device__ void Domain_d::calcDerivatives_FullInt () {
  
  int e = threadIdx.x + blockDim.x*blockIdx.x;
  if (e < m_elem_count) {
    
  // integer :: e
  // ! !rg=gauss[ig]
  // ! !sg=gauss[jg]
  // real(fp_kind), dimension(dim,nodxelem) :: dHrs !!! USED ONLY FOR SEVERAL GAUSS POINTS
  Matrix dHrs; /// IN ELEM_TYPE
  // real(fp_kind), dimension(nodxelem,dim) :: x2
  // real(fp_kind), dimension(dim,dim) :: test
  // real(fp_kind), dimension(dim, dim*nodxelem) :: temph
  
  // integer :: i,j,k, gp
  // real(fp_kind):: r   !!! USED ONLY FOR SEVERAL GAUSS POINTS
  // real(fp_kind), dimension(8,3):: gpc !!! gauss point coordinates, r,s,t
  
  // gp = 1
  // do e=1, elem_count
// ! #ifdef _PRINT_DEBUG_  
    // ! print *, "el ", e 
// ! #endif    
    // do i=1,nodxelem
        // !print *, "elnod " , elem%elnod(e,i)
        // x2(i,:)=nod%x(elem%elnod(e,i),:)
    // end do
    
    // if (elem%gausspc(e) .eq. 1) then      
    
      // if (dim .eq. 2) then 
        // !dHdrs [-1,1,1,-1;  -1.-1,1,1] x X2
        // !! J = [
        // !! dx/dr dy/dr
        // !! dx/ds dy/dx ]
        // !!! THIS IS TO AVOID MATMUL
        // ! print *, "nodes X ", x2(:,1)
        // ! print *, "nodes Y ", x2(:,2)
                
        // elem%jacob(e,gp,1,:) = -x2(1,:)+x2(2,:)+x2(3,:)-x2(4,:)
        // elem%jacob(e,gp,2,:) = -x2(1,:)-x2(2,:)+x2(3,:)+x2(4,:)
        // elem%jacob(e,gp,:,:) = 0.25*elem%jacob(e,gp,:,:)
        // else !!!DIM 3
          // !!!!! SETTING LIKE THIS AVOID MATMUL
          // elem%jacob(e,gp,1,:) = -x2(1,:)+x2(2,:)+x2(3,:)-x2(4,:)-x2(5,:)+x2(6,:)+x2(7,:)-x2(8,:)
          // elem%jacob(e,gp,2,:) = -x2(1,:)-x2(2,:)+x2(3,:)+x2(4,:)-x2(5,:)-x2(6,:)+x2(7,:)+x2(8,:)
          // elem%jacob(e,gp,3,:) = -x2(1,:)-x2(2,:)-x2(3,:)-x2(4,:)+x2(5,:)+x2(6,:)+x2(7,:)+x2(8,:)
          // !elem%jacob(e,gp,2,:) = [-x2(1,2),-x2(2,2), x2(3,2), x2(4,2),-x2(5,2),-x2(6,2), x2(7,2), x2(8,2)]
          // !elem%jacob(e,gp,3,:) = [-x2(1,3),-x2(2,3), x2(3,3), x2(4,3),-x2(5,3),-x2(6,3), x2(7,3), x2(8,3)]
          // ! dHrs(1,:)=[-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0] AND THIS IS dHrs*x2
          // ! dHrs(2,:)=[-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0]       
          // ! dHrs(3,:)=[-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0]  
          // ! elem%jacob(e,gp,1,:) = matmul(dHrs,x2)
          // elem%jacob(e,gp,:,:) = 0.125*elem%jacob(e,gp,:,:)
      // end if  !!!!DIM
      // elem%detJ(e,gp) = det(elem%jacob(e,gp,:,:))
    // else !!!!! GP > 1
      // r = 1.0/sqrt(3.0);
      // gpc(1,:)=[-r,-r,-r];   gpc(2,:)=[ r,-r,-r];      gpc(3,:)=[-r, r,-r];      gpc(4,:)=[ r, r,-r]; !These are the 4 points for 2D full elem
      // gpc(5,:)=[-r,-r, r];   gpc(6,:)=[ r,-r, r];      gpc(7,:)=[-r, r, r];      gpc(8,:)=[ r, r, r];
    
      if (m_dim == 3) {
        for (int gp=0;gp<gp_count;gp++){

          // dHrs(1,:)=[-1.0*(1-gpc(gp,2))*(1.0-gpc(gp,3)),     (1-gpc(gp,2))*(1.0-gpc(gp,3))&
                    // ,     (1+gpc(gp,2))*(1.0-gpc(gp,3)),-1.0*(1+gpc(gp,2))*(1.0-gpc(gp,3))&
                    // ,-1.0*(1-gpc(gp,2))*(1.0+gpc(gp,3)),     (1-gpc(gp,2))*(1.0+gpc(gp,3))&
                    // ,     (1+gpc(gp,2))*(1.0+gpc(gp,3)),-1.0*(1+gpc(gp,2))*(1.0+gpc(gp,3))]
          // dHrs(2,:)=[-1.0*(1-gpc(gp,1))*(1.0-gpc(gp,3)),-1.0*(1+gpc(gp,1))*(1.0-gpc(gp,3))&
                         // ,(1+gpc(gp,1))*(1.0-gpc(gp,3)),     (1-gpc(gp,1))*(1.0-gpc(gp,3))&
                    // ,-1.0*(1-gpc(gp,1))*(1.0+gpc(gp,3)),-1.0*(1+gpc(gp,1))*(1.0+gpc(gp,3))&
                         // ,(1+gpc(gp,1))*(1.0+gpc(gp,3)),     (1-gpc(gp,1))*(1.0+gpc(gp,3))]
          // dHrs(3,:)=[-1.0*(1-gpc(gp,1))*(1.0-gpc(gp,2)),-1.0*(1+gpc(gp,1))*(1.0-gpc(gp,2))&
                    // ,-1.0*(1+gpc(gp,1))*(1.0+gpc(gp,2)),-1.0*(1-gpc(gp,1))*(1.0+gpc(gp,2))&
                    // ,     (1-gpc(gp,1))*(1.0-gpc(gp,2)),     (1+gpc(gp,1))*(1.0-gpc(gp,2))&
                    // ,     (1+gpc(gp,1))*(1.0+gpc(gp,2)),     (1-gpc(gp,1))*(1.0+gpc(gp,2))]                     
          
          // elem%dHrs(e,gp,:,:) =  dHrs(:,:)         
          // !dHrs(2,:)=[(1+r(i)), (1-r(i)),-(1-r(i)),-(1+r(i))]         
          // !dHrs(3,:)=[(1+r(i)), (1-r(i)),-(1-r(i)),-(1+r(i))] 
          // !print *, "dhrs", dHrs 
          // !print *, "x2", x2 
          // elem%jacob(e,gp,:,:) = 0.125*matmul(dHrs,x2)
// ! #if defined _PRINT_DEBUG_
          // ! print *, "jacob ", elem%jacob(e,gp,:,:)
// ! #endif          
          // elem%detJ(e,gp) = det(elem%jacob(e,gp,:,:))
          // !print *, "detJ ", elem%detJ(e,gp)
        }
      } else { //!dim =2
        // do gp = 1,4
          // dHrs(1,:)=[-1.0*(1-gpc(gp,2)),     (1-gpc(gp,2))&
                    // ,     (1+gpc(gp,2)),-1.0*(1+gpc(gp,2))]
          // dHrs(2,:)=[-1.0*(1-gpc(gp,1)),-1.0*(1+gpc(gp,1))&
                         // ,(1+gpc(gp,1)),     (1-gpc(gp,1))]                
          
          // elem%dHrs(e,gp,:,:) =  dHrs(:,:)         
          // !dHrs(2,:)=[(1+r(i)), (1-r(i)),-(1-r(i)),-(1+r(i))]         
          // !dHrs(3,:)=[(1+r(i)), (1-r(i)),-(1-r(i)),-(1+r(i))] 
          // !print *, "dhrs", dHrs 
          // !print *, "x2", x2 
          // elem%jacob(e,gp,:,:) = 0.25*matmul(dHrs,x2)
// ! #if defined _PRINT_DEBUG_
          // !print *, "jacob ", elem%jacob(e,gp,:,:)
// ! #endif          
          // elem%detJ(e,gp) = det(elem%jacob(e,gp,:,:))
          // !print *, "detJ ", elem%detJ(e,gp)
        // end do !gp      
        
      }
    // end if !!gp ==1
// ! #if defined _PRINT_DEBUG_
    // !print *, "jacob ", elem%jacob(e,gp,:,:)
// ! #endif    
  // end do !element
  }
}


};
	