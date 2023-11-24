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

#define __spec __device__ __inline__

class Matrix {
public: 
  __spec Matrix(){}
  __spec Matrix(const int &row, const int &col);
  
  __spec double & getVal(const int &a, const int &b);
  __spec double & operator()(const int &a, const int &b);
  
  __spec ~Matrix(){cudaFree (m_data);}

	double *m_data;
  int m_row, m_col;

};


__spec Matrix::Matrix(const int &row, const int &col) {
  m_row = row;
  m_col = col;
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

#endif
