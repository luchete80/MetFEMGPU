#include "cuda/Domain_d.cuh"

#include <iostream>
#include "cuda/cudautils.cuh"

using namespace MetFEM;

using namespace std;
void report_gpu_mem()
{
    size_t free, total;
    cudaMemGetInfo(&free, &total);
    std::cout << "Free = " << free << " Total = " << total <<std::endl;
}

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

int main(){


	Domain_d *dom_d;

	report_gpu_mem();
	gpuErrchk(cudaMallocManaged(&dom_d, sizeof(MetFEM::Domain_d)) );
	report_gpu_mem();
		
	double3 V = make_double3(0.0,0.0,0.0);
	double3 L = make_double3(0.1,0.1,0.1);
	double r = 0.05;
	
	dom_d->AddBoxLength(V,L,r);
	
	//SolverChungHulbert solver(&dom);
	cout << "Element Count "<<dom_d->getElemCount()<<endl;
	dom_d->SolveChungHulbert ();
	cout << "Program ended."<<endl;
	
	
}