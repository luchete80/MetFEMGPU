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
  
  double *a;
  cudaFree(a);

  ////// MATERIAL  
  double E, nu;
  Material_ *mat_h = (Material_ *)malloc(dom_d->getElemCount() * sizeof(Material_ *)); 
  Elastic_ el(E,nu);
  // cout << "Mat type  "<<mattype<<endl;

  Material_ *material_h;
  double Ep, c[6];
  // MATERIAL
  //TODO: MATERIALS SHOULD BE A VECTOR
  
  string mattype = "Bilinear";
  if      (mattype == "Bilinear")    {
    Ep = E*c[0]/(E-c[0]);		                              //only constant is tangent modulus
    material_h  = new Bilinear(Ep);
    cout << "Material Constants, Et: "<<c[0]<<endl;
    // cudaMalloc((void**)&dom_d->materials, 1 * sizeof(Bilinear )); //
    // cudaMemcpy(dom_d->materials, material_h, 1 * sizeof(Bilinear), cudaMemcpyHostToDevice);	
  } 
  // else if (mattype == "Hollomon")    {
    // // material_h  = new Hollomon(el,Fy,c[0],c[1]);
    // // cout << "Material Constants, K: "<<c[0]<<", n: "<<c[1]<<endl;
    // // cudaMalloc((void**)&dom_d->materials, 1 * sizeof(Hollomon));
    
    // material_h  = new Material_(el);
    // material_h->InitHollomon(el,Fy,c[0],c[1]);
    // material_h->Material_model = HOLLOMON;
    // cudaMalloc((void**)&dom_d->materials, 1 * sizeof(Material_));
    
    // //init_hollomon_mat_kernel<<<1,1>>>(dom_d); //CRASH
    // //cudaMemcpy(dom_d->materials, material_h, 1 * sizeof(Hollomon*), cudaMemcpyHostToDevice);	
    // cudaMemcpy(dom_d->materials, material_h, 1 * sizeof(Material_), cudaMemcpyHostToDevice);	 //OR sizeof(Hollomon)??? i.e. derived class
    
  
  // } else if (mattype == "JohnsonCook") {
    // //Order is 
                               // //A(sy0) ,B,  ,C,   m   ,n   ,eps_0,T_m, T_transition
   // //Material_ *material_h  = new JohnsonCook(el,Fy, c[0],c[1],c[3],c[2],c[6], c[4],c[5]); //First is hardening // A,B,C,m,n_,eps_0,T_m, T_t);	 //FIRST IS n_ than m
    
    // //Only 1 material to begin with
    // //cudaMalloc((void**)&dom_d->materials, 1 * sizeof(JohnsonCook ));
    // //cudaMemcpy(dom_d->materials, material_h, 1 * sizeof(JohnsonCook), cudaMemcpyHostToDevice);	
    // cout << "Material Constants, B: "<<c[0]<<", C: "<<c[1]<<", n: "<<c[2]<<", m: "<<c[3]<<", T_m: "<<c[4]<<", T_t: "<<c[5]<<", eps_0: "<<c[6]<<endl;
  // } else                              printf("ERROR: Invalid material type.

	
	//SolverChungHulbert solver(&dom);
	cout << "Element Count "<<dom_d->getElemCount()<<endl;
	dom_d->SolveChungHulbert ();
	cout << "Program ended."<<endl;
	
	
}