#include "Domain_d.cuh"
#include <iostream>
#include <vector>

#include "tensor.cuh"
// #include "tensor.cu"
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
  Matrix *rot_rate = new Matrix(m_dim,m_dim); //TODO: MAKE SYMM MATRIX
  //Matrix *dHxy_detJ_loc = new Matrix(m_dim, m_nodxelem);

  if (e < m_elem_count) {
    for (int gp=0;gp<m_gp_count;gp++){
      int offset = e * m_gp_count * 6 + gp;
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
        // double vele[3];
        // vele[0] = vele3.x;        vele[1] = vele3.y;        vele[2] = vele3.z;
        // do d=1, dim
          // !print *, "node dim dHxy vele", n,d,temp(d,n) , elem%vele (e,dim*(n-1)+d,1) 
          // elem%str_rate(e,gp, d,d) = elem%str_rate(e,gp, d,d) + temp(d,n) * elem%vele (e,dim*(n-1)+d,1) 
          // elem%rot_rate(e,gp, d,d) = 0.0d0
        // end do
        for (int d=0;d<m_dim;d++){
          str_rate->Set(d,d, str_rate->getVal(d,d) + getDerivative(e,gp,d,n) * f * getVElem(e,n,d));
          // elem%str_rate(e,gp, d,d) = elem%str_rate(e,gp, d,d) + temp(d,n) * elem%vele (e,dim*(n-1)+d,1) 
          // elem%rot_rate(e,gp, d,d) = 0.0d0
        }//dim
        // !!!! TO AVOID ALL MATMULT
        str_rate->Set(0,1, str_rate->getVal(0,1) + f *(getDerivative(e,gp,1,n) * getVElem(e,n,0) +
                                                       getDerivative(e,gp,0,n) * getVElem(e,n,1)));
        rot_rate->Set(0,1, rot_rate->getVal(0,1) + f* (getDerivative(e,gp,1,n) * getVElem(e,n,0) - 
                                                       getDerivative(e,gp,0,n) * getVElem(e,n,1)));
        // elem%str_rate(e,gp, 1,2) = elem%str_rate(e,gp, 1,2) + temp(2,n)* elem%vele (e,dim*(n-1)+1,1) &!!!!dvx/dy
                                   // + temp(1,n) * elem%vele (e,dim*(n-1)+2,1)
        // elem%rot_rate(e,gp, 1,2) = elem%rot_rate(e,gp, 1,2) + temp(2,n)* elem%vele (e,dim*(n-1)+1,1) & !!!!dvx/dx
                                   // - temp(1,n) * elem%vele (e,dim*(n-1)+2,1)                           !!!!
        if (m_dim == 3) {
          str_rate->Set(1,2, str_rate->getVal(1,2) + f *(getDerivative(e,gp,1,n) * getVElem(e,n,0) +
                                                         getDerivative(e,gp,0,n) * getVElem(e,n,1)));
          str_rate->Set(0,2, str_rate->getVal(0,2) + f *(getDerivative(e,gp,2,n) * getVElem(e,n,0) +
                                                         getDerivative(e,gp,0,n) * getVElem(e,n,2)));
          // elem%str_rate(e,gp, 2,3) = elem%str_rate(e,gp, 2,3) + temp(3,n)* elem%vele (e,dim*(n-1)+2,1) &!!!d/dz*vy     
                                     // + temp(2,n) * elem%vele (e,dim*(n-1)+3,1)    !!!d/dy*vz
          // elem%str_rate(e,gp, 1,3) = elem%str_rate(e,gp, 1,3) + temp(3,n)* elem%vele (e,dim*(n-1)+1,1) & !!!d/dz*vx     
                                     // + temp(1,n) * elem%vele (e,dim*(n-1)+3,1)    !!!d/dx*vz     
          // elem%rot_rate(e,gp, 2,3) = elem%rot_rate(e,gp, 2,3) + temp(3,n)* elem%vele (e,dim*(n-1)+2,1) &
                                     // - temp(2,n) * elem%vele (e,dim*(n-1)+3,1)    !!!d/dy*vz
          // elem%rot_rate(e,gp, 1,3) = elem%rot_rate(e,gp, 1,3) + temp(3,n)* elem%vele (e,dim*(n-1)+1,1) & !!!d/dz*vx     
                                     // - temp(1,n) * elem%vele (e,dim*(n-1)+3,1)    !!!d/dx*vz    
        }// end if     
      }// end do !Nod x elem
      
      //str_rate->
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
      str_rate->ToFlatSymPtr(m_str_rate, offset);
      rot_rate->ToFlatSymPtr(m_rot_rate, offset);
      
      printf("Strain Rate\n");
      str_rate->Print();
      
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
    
    delete str_rate,rot_rate;
  } //calcElemStrains
  
  __global__ void calcElemStrainsKernel(Domain_d *dom_d){
		
		dom_d->calcElemStrains();
}

__device__ void Domain_d::calcElemForces(){
  int e = threadIdx.x + blockDim.x*blockIdx.x;  
  if (e < m_elem_count) {
    for (int gp=0;gp<m_gp_count;gp++){
	
  // integer :: e, i,j,k, gp,n, d
  // real(fp_kind), dimension(dim*nodxelem,1) ::f
  // real(fp_kind) :: w
  // !TESTING
  // real (fp_kind) :: sigma_test(6,1) !ORDERED
  // real(fp_kind) :: test(24,1) !ifwanted to test in tensor form
  // elem%f_int = 0.0d0
  // w = 1.0d0 !!! Full integration
	
	// ! !$omp parallel do num_threads(Nproc) private (e,gp,d, w,n) 
  // do e=1, elem_count
    // if (elem%gausspc(e) .eq. 1) then
      // w = 2.0d0**dim
    // end if
    // do gp = 1, elem%gausspc(e)
      // !print *, "elem%dHxy_detJ(e,gp,1", elem%dHxy_detJ(e,gp,1,:)
      // !print *, "elem%dHxy_detJ(e,gp,2", elem%dHxy_detJ(e,gp,2,:)
      // ! sigma_test (:,1)=[elem%sigma (e,gp, 1,1),elem%sigma (e,gp, 2,2),elem%sigma (e,gp, 3,3),&
                        // ! elem%sigma (e,gp, 1,2),elem%sigma (e,gp, 2,3),elem%sigma (e,gp, 3,1)]
      // ! test = w*matmul(transpose(elem%bl(e,gp,:,:)),sigma_test)  ! (24x6)(6x1)
      // !print *, "test force", test
      
      // !print *, "dHdxy, 1", elem%dHxy_detJ(e,gp,1,:)
      // !print *, "dHdxy, 2", elem%dHxy_detJ(e,gp,2,:)
      // !print *, "dHdxy, 3", elem%dHxy_detJ(e,gp,1,:)
      int offset = e*m_nodxelem*m_dim;
      int sig_offset = m_elem_count * m_gp_count * 6;
      
      // do n=1, nodxelem
      // !Is only linear matrix?    
      // !elem%f_int(e,n,d) =  
      // !f (:,:) = matmul(transpose(elem%bl(e,gp,:,:)),elem%sigma (e,:,:))
      // !!!! TO AVOID MATRIX MULTIPLICATIONS (8x6 = 48 in bathe notation with several nonzeros)
      // !!!!! F = BT x sigma = [dh1/dx dh1/dy ] x [ sxx sxy]
      // !!!!!                = [dh2/dx dh2/dy ]   [ syx syy]
      // !!!!! 

      for (int n=0; n<m_nodxelem;n++) {
        for (int d=0;d<m_dim;d++){
          m_f_elem[offset + n*m_dim + d] += getDerivative(e,gp,d,n) * getSigma(e,gp,d,d);
        }
      
        if (m_dim == 2){
          m_f_elem[offset + n*m_dim    ] +=  getDerivative(e,gp,1,n) * getSigma(e,gp,0,1);
          m_f_elem[offset + n*m_dim + 1] +=  getDerivative(e,gp,0,n) * getSigma(e,gp,0,1);
        } else {
          m_f_elem[offset + n*m_dim    ] +=  getDerivative(e,gp,1,n) * getSigma(e,gp,0,1) +
                                             getDerivative(e,gp,2,n) * getSigma(e,gp,0,2);
          m_f_elem[offset + n*m_dim + 1] +=  getDerivative(e,gp,0,n) * getSigma(e,gp,0,1) + 
                                             getDerivative(e,gp,2,n) * getSigma(e,gp,1,2);        
          m_f_elem[offset + n*m_dim + 2] +=  getDerivative(e,gp,1,n) * getSigma(e,gp,1,2) + 
                                             getDerivative(e,gp,0,n) * getSigma(e,gp,0,2);     
        }
      
      }
        // do d=1, dim
          // elem%f_int(e,n,d) = elem%f_int(e,n,d) + elem%dHxy_detJ(e,gp,d,n) * elem%sigma (e,gp, d,d)
        // end do
        // if (dim .eq. 2) then  !!!!! TODO: CHANGE WITH BENSON 1992 - EQ 2.4.2.11 FOR SIMPLICITY
          // !!elem%f_int(e,n,1) = 
          // elem%f_int(e,n,1) = elem%f_int(e,n,1) + elem%dHxy_detJ(e,gp,2,n) * elem%sigma (e,gp, 1,2) 
          // elem%f_int(e,n,2) = elem%f_int(e,n,2) + elem%dHxy_detJ(e,gp,1,n) * elem%sigma (e,gp, 1,2)
        // else 
          // elem%f_int(e,n,1) = elem%f_int(e,n,1) + elem%dHxy_detJ(e,gp,2,n) * elem%sigma (e,gp, 1,2) + &
                                                  // elem%dHxy_detJ(e,gp,3,n) * elem%sigma (e,gp, 1,3)
          // elem%f_int(e,n,2) = elem%f_int(e,n,2) + elem%dHxy_detJ(e,gp,1,n) * elem%sigma (e,gp, 1,2) + &
                                                  // elem%dHxy_detJ(e,gp,3,n) * elem%sigma (e,gp, 2,3)
          // elem%f_int(e,n,3) = elem%f_int(e,n,3) + elem%dHxy_detJ(e,gp,2,n) * elem%sigma (e,gp, 2,3) + &
                                                  // elem%dHxy_detJ(e,gp,1,n) * elem%sigma (e,gp, 1,3)
        // end if
        // !print *, "Element force Node ", n, "F  ", elem%f_int(e,n,:) * w 
      // end do! nod x elem
      // !print *, "test ", w * elem%dHxy_detJ(e,gp,3,8)  * elem%sigma (e,gp, 3,3)
      // !print *, "dHxy ", elem%dHxy_detJ(e,gp,3,8), "w ", w
      // !print *, "s33 ", elem%sigma (e,gp, 3,3)
    // end do !gp
    // elem%f_int(e,:,:) = elem%f_int(e,:,:) * w
  // end do!elem
	// ! !$omp end parallel do    

      } // Gauss Point
    }//if e<elem_count
}

  __global__ void calcElemForcesKernel(Domain_d *dom_d){
		
		dom_d->calcElemForces();
}

__device__ void Domain_d::calcElemPressure(){

  int e = threadIdx.x + blockDim.x*blockIdx.x;  
  if (e < m_elem_count) {
    double press_inc = 0.0;
    for (int gp=0;gp<m_gp_count;gp++){
      for (int d = 0; d<m_dim;d++){
        press_inc += getStrRate(e,gp,d,d); // TODO: TRACE
      }
    }//gauss point
    press_inc = -press_inc/m_gp_count;
// subroutine calc_elem_pressure_from_strain (modK)
  // implicit none
  // real(fp_kind), intent(in) :: modK
  // integer :: e, gp  
  // real(fp_kind) :: press_inc
  
  // gp = 1
  // do e = 1, elem_count
    // press_inc = 0.0d0
    // do gp = 1, elem%gausspc(e)
      // press_inc = press_inc + trace (elem%str_inc(e,gp,:,:))
    // end do
    // press_inc = -press_inc/elem%gausspc(e)
    // ! print *, "press inc ", press_inc
    // do gp = 1, elem%gausspc(e)    
          // elem%pressure(e,gp) = -1.0/3.0 * trace (elem%sigma(e,gp,:,:)) + press_inc * modK
    // end do
    // ! print *, "mod K", modK
    // ! print *, "strain inc", elem%str_inc(e,gp,:,:)
    // ! print *, "press_inc ", press_inc
    // ! print *, "elem%pressure(e,gp) FROM STRAIN", elem%pressure(e,1)
  // end do
// end subroutine 
  } // e< elem_count
}


// //////////////////////////////////////////////////////////////////////
// __device__ void Domain_d::StressStrain(int i) {
	
	// //int i = threadIdx.x + blockDim.x*blockIdx.x;
	// double dep = 0.;
  // double Ep;  //STORAGE THIS OR THE TANGENT?
	
	// if ( i < solid_part_count ) {	
		// //Pressure = EOS(PresEq, Cs, P0,Density, RefDensity); //CALL BEFORE!

		// // Jaumann rate terms
		// tensor3 RotationRateT,SRT,RS;
		// tensor3 RotationRate;
		// tensor3 StrainRate;
		// tensor3 ShearStress;
		// tensor3 Sigma;
		// tensor3 Strain,Straina,Strainb;
		
		// double temprr[6],tempss[6],tempsr[6];
		// double tempssa[6],tempssb[6];
		// for (int k=0;k<6;k++){ //First the diagonal
			// temprr[k] = rotrate[6*i+k];
			// tempss[k] = shearstress[6*i+k];
			// tempsr[k] = strrate[6*i+k];
		// }
		// ShearStress   = FromFlatSym (tempss);	
		// StrainRate    = FromFlatSym(tempsr);
		// RotationRate  = FromFlatAntiSym(temprr);
   // // printf(" Strain rate xx %f \n",StrainRate.zz);
		// RotationRateT = Trans(RotationRate);
		
		// SRT = ShearStress * RotationRateT;
		// RS = RotationRate * ShearStress;

		// ShearStress	= deltat*(2.0*G[i]*(StrainRate-1.0/3.0*(StrainRate.xx+StrainRate.yy+StrainRate.zz)*Identity())+SRT+RS) + ShearStress;	

    // eff_strain_rate[i] = sqrt ( 	0.5*( (StrainRate.xx-StrainRate.yy)*(StrainRate.xx-StrainRate.yy) +
                                // (StrainRate.yy-StrainRate.zz)*(StrainRate.yy-StrainRate.zz) +
                                // (StrainRate.zz-StrainRate.xx)*(StrainRate.zz-StrainRate.xx)) + 
                          // 3.0 * (StrainRate.xy*StrainRate.xy + StrainRate.yz*StrainRate.yz + StrainRate.zx*StrainRate.zx)
                        // );
                        
		// double J2	= 0.5*(ShearStress.xx*ShearStress.xx + 2.0*ShearStress.xy*ShearStress.yx +
					// 2.0*ShearStress.xz*ShearStress.zx + ShearStress.yy*ShearStress.yy +
					// 2.0*ShearStress.yz*ShearStress.zy + ShearStress.zz*ShearStress.zz);

    // //Scale back, Fraser Eqn 3-53
		// double sig_trial = sqrt(3.0*J2); 
    // if ( sigma_y[i] < sig_trial ) ShearStress = sigma_y[i]/sig_trial * ShearStress; //Yielding      
    // //std::min((Sigmay/sqrt(3.0*J2)),1.0)*ShearStressa;

   // if       (mat[i]->Material_model == HOLLOMON )       {
     
      // //printf("Hollomon!");
      // //printf("pl strain %f\n",pl_strain[i]);
      // //sigma_y[i] = mat [i]->CalcYieldStress(pl_strain[i]);
      // //ShowProps(mat[i]);
      // sigma_y[i] = CalcHollomonYieldStress(pl_strain[i],mat [i]);
      // //printf ("sigma_y %f\n",sigma_y[i]);
      // //sigma_y[i] = mat [i]->testret();
      // //(*materials_ptr)->testret();
      // //sigma_y[i] = mat [i]->CalcYieldStress(pl_strain[i]); 
    // } 
		// else if  (mat[i]->Material_model == JOHNSON_COOK )   sigma_y[i] = CalcJohnsonCookYieldStress(pl_strain[i],eff_strain_rate[i],T[i], mat[i]);
		
		// sigma_eq[i] = sig_trial;	
		
		// if ( sig_trial > sigma_y[i]) {
      // if              (mat[i]->Material_model == HOLLOMON ) {
				// Et[i] = CalcHollomonTangentModulus(pl_strain[i],mat[i]); //Fraser 3.54
				// // Et_m = Et;        
      // } else if       (mat[i]->Material_model == JOHNSON_COOK ) {
        // Et[i] = CalcJohnsonCookTangentModulus(pl_strain[i], eff_strain_rate[i], T[i],mat[i]); //Fraser 3.54
      // } else if       (mat[i]->Material_model == BILINEAR ) {
        // //Ep = mat->Elastic().E()*Et/(mat->Elastic().E()-Et);
      // }
			// if (mat[i]->Material_model > BILINEAR ) {//Else Ep = 0
        // //cout << "Calculating Ep"<<endl;
				// Ep = mat[i]->Elastic().E()*Et[i]/(mat[i]->Elastic().E()-Et[i]);
				// // if (Ep < 0)
					// // cout << "ATTENTION Material Ep <0 "<<Ep<<", Et" << Et <<", platrain"<<pl_strain<<"effstrrate"<<eff_strain_rate<<endl;
			// }
      // if (Ep<0) Ep = 1.*mat[i]->Elastic().E();
      
			// dep=( sig_trial - sigma_y[i])/ (3.*G[i] /*+ Ep*/);	//Fraser, Eq 3-49 TODO: MODIFY FOR TANGENT MODULUS = 0
			// pl_strain[i] += dep;	
      // //printf("Particle %d, dep %.1e, sigtrial %.1e\n",i,dep,sig_trial);
			// sigma_eq[i] = sigma_y[i];
		// }

    // // if (mat[i]->Material_model == BILINEAR )
      // // sigma_y[i] += dep*mat[i]->Ep;
    
		// Sigma = -p[i] * Identity() + ShearStress;	//Fraser, eq 3.32

		// Strain	= deltat * StrainRate + Strain;
    
    // // if (mat[i]->Material_model==JOHNSON_COOK){
      // // printf("JOHNSON_COOK!\n"); //test
    // // } elsr     
    // // if (mat[i]->Material_model==HOLLOMON){
      // // printf("HOLLOMON!\n"); //test
    // // }    
    // // if (mat[i]->Material_model==BILINEAR){
      // // printf("BILINEAR!\n"); //test
    // // }

		// ///// OUTPUT TO Flatten arrays
		// ToFlatSymPtr(Sigma, sigma,6*i);  //TODO: CHECK IF RETURN VALUE IS SLOWER THAN PASS AS PARAM		
		// ToFlatSymPtr(Strain, 	strain,6*i);		
		// ToFlatSymPtr(ShearStress, shearstress,6*i);

	// }//particle count
// }

  __global__ void calcElemPressureKernel(Domain_d *dom_d){
		
		dom_d->calcElemPressure();
}

 
}; //Namespace MetFEM