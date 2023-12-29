#include "Domain_d.cuh"

#include <cuda_runtime.h>

namespace MetFEM {
  
// __global__ void histogram(int* color, int* buckets)
// {
 // int i = threadIdx.x + blockDim.x * blockIdx.x;
 // int c = colors[i];
 // atomicAdd(&buckets[c], 1);
// }

__device__ void Domain_d::assemblyForces(){

  int e = threadIdx.x + blockDim.x*blockIdx.x;
  if (e < m_elem_count) {
    for (int n=0; n<m_nodxelem;n++) {
      for (int d=0;d<m_dim;d++){
        atomicAdd(&m_f[m_elnod[n]*m_dim + d], m_f_elem[e*m_nodxelem*m_dim + n*m_dim + d]);
      }
    }
  } // element
  
// subroutine assemble_forces()
  // integer :: e, i, n, iglob
  // real(fp_kind), dimension(nodxelem*dim,1) :: utemp, rtemp
  
  // !print *, "assemblying int forces"
  // rint_glob (:,:) = 0.0d0
  // fext_glob (:,:) = 0.0d0
  // do e = 1, elem_count
    // do n = 1, nodxelem
      // !print *,"elem mat kl", elem%matkl(e,:,:)
      // !print *, "elem fext ", elem%f_ext(e,n,:)
      // !print *, "elem fint ", elem%f_int(e,n,:)
      // do i=1,dim 
        // iglob  = dim * (elem%elnod(e,n) - 1 ) + i
        // !rint_glob(elem%elnod(e,n),i) =  rint_glob(elem%elnod(e,n),i) + rtemp(dim*(n-1)+i,1)
        // !print *, "rint ",rint_glob(elem%elnod(e,n),i), "fint ", elem%f_int(e,n,i)
        // rint_glob(elem%elnod(e,n),i) =  rint_glob(elem%elnod(e,n),i) + elem%f_int(e,n,i)
        // fext_glob(elem%elnod(e,n),i) =  fext_glob(elem%elnod(e,n),i) + elem%f_ext(e,n,i)
      // end do !element row
    // end do ! Element node
    // if (elem%gausspc(e) .eq. 1) then
      // do n = 1, nodxelem
        // rint_glob(elem%elnod(e,n),:) = rint_glob(elem%elnod(e,n),:) + elem%hourg_nodf(e,n,:)
      // end do
    // end if 
  // end do ! e
// end subroutine 

}

}