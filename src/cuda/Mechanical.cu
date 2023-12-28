#include "Domain_d.cuh"
#include <iostream>
#include <vector>

#include "tensor.cuh"
#include "Matrix.cuh"

using namespace std;

namespace MetFEM {
  
  
// !!!!!!!!!!!!!!!Gradv = L = dvx/dx dvx/dy  dvx/dz
// !!!!!!!!!!!!!!!!!!!        dvy/dx dvy/dy  dvy/dz
// !!!! E = 1/2 (L+LT)
// !!!! R = 1/2 (L-LT)
// !THIS SHOULD BE DONE AT t+1/2dt
// subroutine cal_elem_strains ()
  // implicit none
  // integer :: e, i,j,k, gp, d, n
  // real(fp_kind), dimension(dim,nodxelem) ::temp
  // real(fp_kind) :: f
  // real(fp_kind) :: test(1,6),test33(3,3) !ifwanted to test in tensor form
  
  // elem%str_rate = 0.0d0
  // elem%rot_rate = 0.0d0
  
  // do e=1, elem_count
    // do gp = 1, elem%gausspc(e)
      // !Is only linear matrix?    
      // !!!TODO: CHANGE FROM MATRIX OPERATION TO SIMPLE OPERATION
      // f = 1.0d0/elem%detJ(e,gp)
      // temp = elem%dHxy_detJ(e,gp,:,:) * f!!!!TODO: MODIFY BY MULTIPLYING
      // elem%strain(e,gp,:,:) = matmul(elem%bl(e,gp,:,:),elem%uele (e,:,:)) 
      // !print *, "standard stran rate calc (matricial) "
      // ! !!!! DEFAULT (TODO: CHECK IF IS SLOW)
      // test = f* matmul(elem%bl(e,gp,:,:),elem%vele (e,:,:))  ! (6x24)(24x1)
      // !print *, "e11 e22 e33 2e12 2e23 2e31", test
      // test33(1,1) = test(1,1);test33(2,2) = test(1,2);test33(3,3) = test(1,3);
      // test33(1,2) = test(1,4)*0.5;test33(2,1) =test33(1,2);
      // test33(2,3) = test(1,5)*0.5;test33(3,2) =test33(2,3);
      // test33(3,1) = test(1,6)*0.5;test33(1,3) =test33(3,1);
      

      // test33 = 0.5*(test33+transpose(test33));
      // !print *, "str rate test", test33
      
      // ! test33 = 0.5*(test33-transpose(test33));
      // ! print *, "rot rate test", test33

      // do n=1, nodxelem  
        // do d=1, dim
          // !print *, "node dim dHxy vele", n,d,temp(d,n) , elem%vele (e,dim*(n-1)+d,1) 
          // elem%str_rate(e,gp, d,d) = elem%str_rate(e,gp, d,d) + temp(d,n) * elem%vele (e,dim*(n-1)+d,1) 
          // elem%rot_rate(e,gp, d,d) = 0.0d0
        // end do
        // !!!! TO AVOID ALL MATMULT
        // elem%str_rate(e,gp, 1,2) = elem%str_rate(e,gp, 1,2) + temp(2,n)* elem%vele (e,dim*(n-1)+1,1) &!!!!dvx/dy
                                   // + temp(1,n) * elem%vele (e,dim*(n-1)+2,1)
        // elem%rot_rate(e,gp, 1,2) = elem%rot_rate(e,gp, 1,2) + temp(2,n)* elem%vele (e,dim*(n-1)+1,1) & !!!!dvx/dx
                                   // - temp(1,n) * elem%vele (e,dim*(n-1)+2,1)                           !!!!
        // if (dim == 3) then
          // elem%str_rate(e,gp, 2,3) = elem%str_rate(e,gp, 2,3) + temp(3,n)* elem%vele (e,dim*(n-1)+2,1) &!!!d/dz*vy     
                                     // + temp(2,n) * elem%vele (e,dim*(n-1)+3,1)    !!!d/dy*vz
          // elem%str_rate(e,gp, 1,3) = elem%str_rate(e,gp, 1,3) + temp(3,n)* elem%vele (e,dim*(n-1)+1,1) & !!!d/dz*vx     
                                     // + temp(1,n) * elem%vele (e,dim*(n-1)+3,1)    !!!d/dx*vz     
          // elem%rot_rate(e,gp, 2,3) = elem%rot_rate(e,gp, 2,3) + temp(3,n)* elem%vele (e,dim*(n-1)+2,1) &
                                     // - temp(2,n) * elem%vele (e,dim*(n-1)+3,1)    !!!d/dy*vz
          // elem%rot_rate(e,gp, 1,3) = elem%rot_rate(e,gp, 1,3) + temp(3,n)* elem%vele (e,dim*(n-1)+1,1) & !!!d/dz*vx     
                                     // - temp(1,n) * elem%vele (e,dim*(n-1)+3,1)    !!!d/dx*vz    
        // end if     
      // end do !Nod x elem
      // elem%str_rate(e,gp, 1,2) = 0.5 * elem%str_rate(e,gp, 1,2); 
      // elem%rot_rate(e,gp, 1,2) = 0.5 * elem%rot_rate(e,gp, 1,2)      

      // elem%str_rate(e,gp, 2,1) =     elem%str_rate(e,gp, 1,2)
      // elem%rot_rate(e,gp, 2,1) =    -elem%rot_rate(e,gp, 1,2)
      // if (dim .eq. 3) then
        // elem%str_rate(e,gp, 1,3) = 0.5 * elem%str_rate(e,gp, 1,3); elem%str_rate(e,gp, 2,3) = 0.5 * elem%str_rate(e,gp, 2,3)
        // elem%rot_rate(e,gp, 1,3) = 0.5 * elem%rot_rate(e,gp, 1,3); elem%rot_rate(e,gp, 2,3) = 0.5 * elem%rot_rate(e,gp, 2,3)
        
        // elem%str_rate(e,gp, 3,2) =     elem%str_rate(e,gp, 2,3)
        // elem%str_rate(e,gp, 3,1) =     elem%str_rate(e,gp, 1,3)

        // elem%rot_rate(e,gp, 3,2) =     -elem%rot_rate(e,gp, 2,3)
        // elem%rot_rate(e,gp, 3,1) =     -elem%rot_rate(e,gp, 1,3)
      // end if

      // !elem%str_rate(e,gp,:,:) = matmul(elem%bl(e,gp,:,:),elem%vele (e,:,:)) 
      // !print *, "simlpified strain rate "
      // !print *, "strain rate ", elem%str_rate(e,gp,:,:)
      // !print *, "rot    rate ", elem%rot_rate(e,gp,:,:)
    // end do !gp
  // end do !element
// end subroutine

  __device__ void Domain_d::calcElemStrains(){

  int e = threadIdx.x + blockDim.x*blockIdx.x;
  Matrix *str_rate = new Matrix(m_dim,m_dim); //TODO: MAKE SYMM MATRIX
  //Matrix *dHxy_detJ_loc = new Matrix(m_dim, m_nodxelem);

  if (e < m_elem_count) {
    for (int gp=0;gp<m_gp_count;gp++){
      int offset = e * m_gp_count;
  // elem%str_rate = 0.0d0
  // elem%rot_rate = 0.0d0
  
      // for (int i=0;i<m_nodxelem;i++ ) {
        // dHxy_detJ_loc->Set(0,i,m_dH_detJ_dx[e*offset + gp * m_gp_count + i]);
        // dHxy_detJ_loc->Set(1,i,m_dH_detJ_dx[e*offset + gp * m_gp_count + i]);
        // dHxy_detJ_loc->Set(2,i,m_dH_detJ_dx[e*offset + gp * m_gp_count + i]);
      // }
      // f = 1.0d0/elem%detJ(e,gp)
      // temp = elem%dHxy_detJ(e,gp,:,:) * f!!!!TODO: MODIFY BY MULTIPLYING
      // elem%strain(e,gp,:,:) = matmul(elem%bl(e,gp,:,:),elem%uele (e,:,:)) 
      // !print *, "standard stran rate calc (matricial) "

      double f = 1.0 / m_detJ[offset + gp];
      for (int n=0; n<m_nodxelem;n++) {
        double3 vele3 = getVElem(e,n);
        double vele[3];
        vele[0] = vele3.x;        vele[1] = vele3.y;        vele[2] = vele3.z;
        // do d=1, dim
          // !print *, "node dim dHxy vele", n,d,temp(d,n) , elem%vele (e,dim*(n-1)+d,1) 
          // elem%str_rate(e,gp, d,d) = elem%str_rate(e,gp, d,d) + temp(d,n) * elem%vele (e,dim*(n-1)+d,1) 
          // elem%rot_rate(e,gp, d,d) = 0.0d0
        // end do
        for (int d=0;d<m_dim;d++){
          str_rate->Set(d,d, str_rate->getVal(d,d) + getDerivative(e,gp,d,n) * f * vele[d]);
          // elem%str_rate(e,gp, d,d) = elem%str_rate(e,gp, d,d) + temp(d,n) * elem%vele (e,dim*(n-1)+d,1) 
          // elem%rot_rate(e,gp, d,d) = 0.0d0
        }//dim
        // !!!! TO AVOID ALL MATMULT
        // elem%str_rate(e,gp, 1,2) = elem%str_rate(e,gp, 1,2) + temp(2,n)* elem%vele (e,dim*(n-1)+1,1) &!!!!dvx/dy
                                   // + temp(1,n) * elem%vele (e,dim*(n-1)+2,1)
        // elem%rot_rate(e,gp, 1,2) = elem%rot_rate(e,gp, 1,2) + temp(2,n)* elem%vele (e,dim*(n-1)+1,1) & !!!!dvx/dx
                                   // - temp(1,n) * elem%vele (e,dim*(n-1)+2,1)                           !!!!
        // if (dim == 3) then
          // elem%str_rate(e,gp, 2,3) = elem%str_rate(e,gp, 2,3) + temp(3,n)* elem%vele (e,dim*(n-1)+2,1) &!!!d/dz*vy     
                                     // + temp(2,n) * elem%vele (e,dim*(n-1)+3,1)    !!!d/dy*vz
          // elem%str_rate(e,gp, 1,3) = elem%str_rate(e,gp, 1,3) + temp(3,n)* elem%vele (e,dim*(n-1)+1,1) & !!!d/dz*vx     
                                     // + temp(1,n) * elem%vele (e,dim*(n-1)+3,1)    !!!d/dx*vz     
          // elem%rot_rate(e,gp, 2,3) = elem%rot_rate(e,gp, 2,3) + temp(3,n)* elem%vele (e,dim*(n-1)+2,1) &
                                     // - temp(2,n) * elem%vele (e,dim*(n-1)+3,1)    !!!d/dy*vz
          // elem%rot_rate(e,gp, 1,3) = elem%rot_rate(e,gp, 1,3) + temp(3,n)* elem%vele (e,dim*(n-1)+1,1) & !!!d/dz*vx     
                                     // - temp(1,n) * elem%vele (e,dim*(n-1)+3,1)    !!!d/dx*vz    
        // end if     
      }// end do !Nod x elem
      
      // elem%str_rate(e,gp, 1,2) = 0.5 * elem%str_rate(e,gp, 1,2); 
      // elem%rot_rate(e,gp, 1,2) = 0.5 * elem%rot_rate(e,gp, 1,2)      

      // elem%str_rate(e,gp, 2,1) =     elem%str_rate(e,gp, 1,2)
      // elem%rot_rate(e,gp, 2,1) =    -elem%rot_rate(e,gp, 1,2)
      // if (dim .eq. 3) then
        // elem%str_rate(e,gp, 1,3) = 0.5 * elem%str_rate(e,gp, 1,3); elem%str_rate(e,gp, 2,3) = 0.5 * elem%str_rate(e,gp, 2,3)
        // elem%rot_rate(e,gp, 1,3) = 0.5 * elem%rot_rate(e,gp, 1,3); elem%rot_rate(e,gp, 2,3) = 0.5 * elem%rot_rate(e,gp, 2,3)
        
        // elem%str_rate(e,gp, 3,2) =     elem%str_rate(e,gp, 2,3)
        // elem%str_rate(e,gp, 3,1) =     elem%str_rate(e,gp, 1,3)

        // elem%rot_rate(e,gp, 3,2) =     -elem%rot_rate(e,gp, 2,3)
        // elem%rot_rate(e,gp, 3,1) =     -elem%rot_rate(e,gp, 1,3)
      // end if

      // !elem%str_rate(e,gp,:,:) = matmul(elem%bl(e,gp,:,:),elem%vele (e,:,:)) 
      // !print *, "simlpified strain rate "

      //Inverse test
      // Matrix *test = new Matrix(3,3);
      // Matrix *invtest = new Matrix(3,3);
      // printf("A\n");
      // test->Set(0,0,1);test->Set(0,1,1);test->Set(0,2,1);
      // test->Set(1,0,1);test->Set(1,1,2);test->Set(1,2,2);      
      // test->Set(2,0,1);test->Set(2,1,2);test->Set(2,2,3);
      // InvMat(*test,invtest);
      // //printf("inv A\n");
      // test->Print();
      // invtest->Print();
        // delete test, invtest;
      } // Gauss Point
    }//if e<elem_count
    
    delete str_rate;
  } //calcElemStrains
  
  __global__ void calcElemStrainsKernel(Domain_d *dom_d){
		
		dom_d->calcElemStrains();
}

  
  
}; //Namespace MetFEM