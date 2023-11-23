#include "cuda/Domain_d.cuh"

using namespace MetFEM;

int main(){

	Domain_d dom;
	double3 V = make_double3(0.0,0.0,0.0);
	double3 L = make_double3(0.1,0.012,0.006);
	double r = 0.0005;
	
	dom.AddBoxLength(V,L,r);
	
}