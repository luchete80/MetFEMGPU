#include "cuda/Domain_d.cuh"

#include <iostream>

using namespace MetFEM;

using namespace std;

int main(){

	Domain_d dom;
  
	double3 V = make_double3(0.0,0.0,0.0);
	double3 L = make_double3(0.1,0.1,0.1);
	double r = 0.02;
	
	dom.AddBoxLength(V,L,r);
	
	//SolverChungHulbert solver(&dom);
	
	dom.SolveChungHulbert ();
	cout << "Program ended."<<endl;
	
	
}