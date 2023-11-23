#include "Domain_d.cuh"
#include <iostream>
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
    
    nel[0] = (int)(L.x/(2.0*r));
    nel[1] = (int)(L.y/(2.0*r));
    cout << "Nel x: "<<nel[0]<<", y "<<nel[1]<<endl;
    if (m_dim == 2){
      nel[2] = 0;
      nodxelem = 4;
    } else {
      nel[2] = (int)(L.z/(2.0*r));
      nodxelem = 8;
    }
    

    Xp.z = V.z ;
    

    // write (*,*) "Creating Mesh ...", "Elements ", neL.y, ", ",neL.z
  int nc = (nel[0] +1) * (nel[1]+1) * (nel[2]+1);
  int ne = nel[0]*nel[1]*nel[2];
  //thisAllocateNodes((neL.x +1) * (nel[1]+1) * (nel[2]+1));
    // print *, "Element count in XYZ: ", nel(:)
    // write (*,*) "Box Node count ", node_count

	this->SetDimension(nc,ne);	 //AFTER CREATING DOMAIN
  //SPH::Domain	dom;
	//double3 *x =  (double3 *)malloc(dom.Particles.size());
	double3 *x =  new double3 [m_node_count];
	for (int i=0;i<dom.Particles.size();i++){
		//cout <<"i; "<<i<<endl;
		//x[i] = make_double3(dom.Particles[i]->x);
		x[i] = make_double3(double(dom.Particles[i]->x(0)), double(dom.Particles[i]->x(1)), double(dom.Particles[i]->x(2)));
	}
	int size = dom.Particles.size() * sizeof(double3);
	cout << "Copying to device..."<<endl;
	cudaMemcpy(this->x, x, size, cudaMemcpyHostToDevice);    
    
    // write (*,*) "xp ", Xp(:)    
    
    if (m_dim == 2) {
    cout << "Box Particle Count is " << m_node_count <<endl;
    p = 1;
      Xp.y = V(1);
      for (int j = 0; j < (nel[1] +1);j++){
        Xp.x = V(0);
        for (int i = 0; i < (neL.x +1);i++){
					m_node.push_back(new Node(Xp));
          //nod%x(p,:) = Xp(:);
          cout << "node " << p <<"X: "<<Xp<<endl;
          p++;
          Xp.x = Xp.x + 2 * r;
        }
        Xp.y = Xp.y + 2 * r;
      }// 
      Xp(2) = Xp(2) + 2 * r;

    cout <<"m_node size"<<m_node.size()<<endl;
    } else {
      // p = 1
      // k = 1; Xp(3) = V(3)
      // do while (k <= (nel(3) +1))
        // j = 1;         Xp(2) = V.z
        // do while (j <= (neL.z +1))
          // i = 1
          // Xp.y = V(1)
          // do while (i <= (neL.y +1))
            // nod%x(p,:) = Xp(:)
            // print *,"node ",p , "X: ",Xp(:)
            // p = p + 1
            // Xp.y = Xp.y + 2 * r
            // i = i +1
          // end do
          // Xp(2) = Xp(2) + 2 * r
          // j = j +1
        // end do 
        // Xp(3) = Xp(3) + 2 * r
        // k = k + 1
      // end do    
		}//if dim
    
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
    
		int ex, ey, ez;
		
    if (m_dim == 2) {
			n.resize(4);
      for (int ey = 0; ey < nel[1];ey++){

           for (int ex = 0; ex < neL.x;ex++){
        int iv[4];
        iv[0] = (neL.x+1)*ey + ex;        iv[1] = (neL.x+1)*ey + ex+1;
        iv[2] = (neL.x+1)*(ey+1) + ex+1;        iv[3] = (neL.x+1)*(ey+1) + ex;
        // cout << i[]
						n[0]= m_node[iv[0]];
						n[1]= m_node[(neL.x+1)*ey + ex+1];
						n[2]= m_node[(neL.x+1)*(ey+1)+ex+1];
						n[3]= m_node[(neL.x+1)*(ey+1)+ex];
            cout << "Nel x : "<<neL.x<<endl;
           cout << "nodes "<<endl;
           for (int i=0;i<4;i++)cout << iv[i]<<", ";
						 m_element.push_back(new El4N2DPE(n));
																							// m_node[(neL.x+1)*ey + ex+1],
																							// m_node[(neL.x+1)*(ey+1)+ex+1],
																							// m_node[(neL.x+1)*(ey+1)+ex]
																							// );
              //elem%elnod(i,:)=[(neL.y+1)*ey + ex+1,(neL.y+1)*ey + ex+2,(neL.y+1)*(ey+1)+ex+2,(neL.y+1)*(ey+1)+ex+1]         
              //print *, "Element ", i , "Elnod", elem%elnod(i,:) 
					 }
      } 
    } else {
      // ez = 0
      // i = 1
      // nnodz = (neL.y+1)*(neL.z+1)
      // print *, "Element Nodes at z ", nnodz
      // do while ( ez < nel(3))
        // ey = 0    
        // do while ( ey < neL.z)
            // ex = 0
            // do while (ex < neL.y) 
                // !elem%elnod(i,:)=[(neL.y+1)*(ey+1)+ex+2,(neL.y+1)*(ey+1)+ex+1,(neL.y+1)*ey + ex+1,(neL.y+1)*ey + ex+2]       
                // elem%elnod(i,:) = [ nnodz*ez + (neL.y+1)*ey + ex+1,nnodz*ez + (neL.y+1)*ey + ex+2, &
                                    // nnodz*ez + (neL.y+1)*(ey+1)+ex+2,nnodz*ez + (neL.y+1)*(ey+1)+ex+1, &
                                    // nnodz*(ez + 1) + (neL.y+1)*ey + ex+1,nnodz*(ez + 1) + (neL.y+1)*ey + ex+2, &
                                    // nnodz*(ez + 1) + (neL.y+1)*(ey+1)+ex+2,nnodz*(ez + 1)+ (neL.y+1)*(ey+1)+ex+1]
                // print *, "Element ", i , "Elnod", elem%elnod(i,:) 
                // i=i+1
              // ex = ex + 1
            // end do
          // ey = ey + 1
        // end do 
        // ez = ez + 1
      // end do !el z
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
}


};
	