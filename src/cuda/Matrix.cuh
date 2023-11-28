/*  Copyright (c) 2013-2015 INGV, EDF, UniCT, JHU

    Istituto Nazionale di Geofisica e Vulcanologia, Sezione di Catania, Italy
    Électricité de France, Paris, France
    Università di Catania, Catania, Italy
    Johns Hopkins University, Baltimore (MD), USA

    This file is part of GPUSPH. Project founders:
        Alexis Hérault, Giuseppe Bilotta, Robert A. Dalrymple,
        Eugenio Rustico, Ciro Del Negro
    For a full list of authors and project partners, consult the logs
    and the project website <https://www.gpusph.org>

    GPUSPH is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GPUSPH is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GPUSPH.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _MATRIX_CUH_
#define _MATRIX_CUH_

#include <cuda.h>
#define __spec __device__ __inline__

class Matrix {
public: 
  __spec Matrix(){}
  __spec Matrix(const int &row, const int &col);
  
  __spec double & getVal(const int &a, const int &b);
  __spec double & operator()(const int &a, const int &b);
	
	__spec void Print();
  
  __spec ~Matrix(){/*cudaFree (m_data);*/}
	
	__spec double calcDet();

	double *m_data;
  int m_row, m_col, m_dim;

};


__spec Matrix::Matrix(const int &row, const int &col) {
  m_row = row;
  m_col = col;
	if (m_row == m_col) m_dim = m_row;
  cudaMalloc((void**)&m_data, row * col * sizeof(double));
  for (int i=0;i<row*col;i++) m_data[i] = 0.0;
}

__spec Matrix MatMul(Matrix &A, Matrix &B){
  Matrix ret(A.m_row,B.m_col);
  for (int i = 0; i<A.m_row; i++)
    for (int j = 0; j<A.m_col; j++)
      for (int k = 0; k<A.m_col; k++)
        ret.m_data[i * A.m_row + j] += A.m_data[i * A.m_row + k] * B.m_data[k * B.m_row + j ];
  
  
  
  return ret;
}

  __spec double & Matrix::getVal(const int &a, const int &b){
    return m_data[m_row*a+b];
  }
  
  __spec double & Matrix::operator()(const int &a, const int &b){
    return m_data[m_row*a+b];
  }
	
	__spec Matrix operator*(const double &c, Matrix &A) {
	Matrix ret;
  for (int i=0;i<A.m_row*A.m_col;i++) ret.m_data[i] = A.m_data[i] * c;
	return ret;
}

	__spec void Matrix::Print() {

  for (int i=0;i<m_row*m_col;i++) {
		for (int j=0;j<m_col*m_col;j++) 
			printf("%lf ", getVal(i,j)) ;
		printf("\n");
	}

}



// subroutine M33INV (A, AINV, OK_FLAG)

  // IMPLICIT NONE

  // DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
  // DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV
  // LOGICAL, INTENT(OUT) :: OK_FLAG

  // DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
  // DOUBLE PRECISION :: DETE
  // DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR


  // ! DET =   A(1,1)*A(2,2)*A(3,3)  &
        // ! - A(1,1)*A(2,3)*A(3,2)  &
        // ! - A(1,2)*A(2,1)*A(3,3)  &
        // ! + A(1,2)*A(2,3)*A(3,1)  &
        // ! + A(1,3)*A(2,1)*A(3,2)  &
        // ! - A(1,3)*A(2,2)*A(3,1)
  
  // dete = det(A)

  // IF (ABS(DETE) .LE. EPS) THEN
     // AINV = 0.0D0
     // OK_FLAG = .FALSE.
     // RETURN
  // END IF

  // COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
  // COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
  // COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
  // COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
  // COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
  // COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
  // COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
  // COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
  // COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

  // AINV = TRANSPOSE(COFACTOR) / DETE

  // OK_FLAG = .TRUE.

// RETURN

// end subroutine M33INV

__spec double Matrix::calcDet (){
	double ret = 0.0;
  // real(fp_kind), dimension(dim,dim), intent (in) :: a 
  // real(fp_kind) :: det
  // if (dim .eq. 2) then
    // det = a(1,1)*a(2,2)-a(1,2)*a(2,1)
  // else 
  // DET =   A(1,1)*A(2,2)*A(3,3)  &
        // - A(1,1)*A(2,3)*A(3,2)  &
        // - A(1,2)*A(2,1)*A(3,3)  &
        // + A(1,2)*A(2,3)*A(3,1)  &
        // + A(1,3)*A(2,1)*A(3,2)  &
        // - A(1,3)*A(2,2)*A(3,1)  
  // end if
	if (m_dim == 2) {
		ret = getVal(0,0) * getVal(1,1) - getVal(0,1) * getVal(1,0);
	} else if (m_dim ==3) {
		ret =   getVal(0,0) * getVal(1,1) * getVal(2,2)
          - getVal(0,0) * getVal(1,2) * getVal(2,1)
					- getVal(0,1) * getVal(1,0) * getVal(2,2)
					+ getVal(0,1) * getVal(1,2) * getVal(2,0)
          + getVal(0,2) * getVal(1,0) * getVal(2,1)
					- getVal(0,2) * getVal(1,1) * getVal(2,0);


	}
	return ret;
}

// function invmat (a)
  // real(fp_kind), dimension(dim,dim), intent (in) :: a 
  // real(fp_kind), dimension(dim,dim) :: invmat 
  // !if (dim .eq. 2) then
  // invmat(1,:) = 1.0d0/(det(a))*[ a(2,2),-a(1,2)]
  // invmat(2,:) = 1.0d0/(det(a))*[-a(2,1), a(1,1)]
  // !end if
// end function

// function adj (a)
  // real(fp_kind), dimension(dim,dim), intent (in) :: a 
  // real(fp_kind), dimension(dim,dim) :: cofactor,adj
  
  // if (dim .eq. 2) then
    // adj(1,:) = [ a(2,2),-a(1,2)]
    // adj(2,:) = [-a(2,1), a(1,1)]
  // else
    // cofactor(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
    // cofactor(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
    // cofactor(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
    // cofactor(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
    // cofactor(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
    // cofactor(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
    // cofactor(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
    // cofactor(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
    // cofactor(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))
    
    // adj = TRANSPOSE(COFACTOR)
  // end if
// end function

#endif
