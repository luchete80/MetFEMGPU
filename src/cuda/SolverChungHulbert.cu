#include "Solver.cuh"
#include <iostream>

using namespace std;

namespace MetFEM{


SolverChungHulbert::SolverChungHulbert(Domain_d *d){
	m_dom = d;
}

__host__ void SolverChungHulbert::Solve(){
	
	int N = m_dom->getElemCount();
	m_dom->threadsPerBlock = 256; //Or BlockSize
	//m_dom->threadsPerBlock = 1; //Or BlockSize
	m_dom->blocksPerGrid =				// Or gridsize
	(N + m_dom->threadsPerBlock - 1) / m_dom->threadsPerBlock;
	cout << "Blocks per grid"<<m_dom->blocksPerGrid<<", Threads per block"<< m_dom->threadsPerBlock<<endl;
	
	calcElemJAndDerivKernel<<<m_dom->blocksPerGrid,m_dom->threadsPerBlock >>>(m_dom);
	cudaDeviceSynchronize(); 
  
  calcElemStrains<<<m_dom->blocksPerGrid,m_dom->threadsPerBlock >>>(m_dom);
	cudaDeviceSynchronize(); 
	
}

}; //Namespace