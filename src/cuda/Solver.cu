#include "Domain_d.cuh"
#include <iostream>

using namespace std;

namespace MetFEM{

	void __host__ Domain_d::SolveChungHulbert(){
	
	int N = getElemCount();
	threadsPerBlock = 256; //Or BlockSize
	//threadsPerBlock = 1; //Or BlockSize
	blocksPerGrid =				// Or gridsize
	(N + threadsPerBlock - 1) / threadsPerBlock;
	cout << "Blocks per grid"<<blocksPerGrid<<", Threads per block"<< threadsPerBlock<<endl;
	
	calcElemJAndDerivKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	cudaDeviceSynchronize(); 
	}
	
};