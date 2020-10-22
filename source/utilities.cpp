#include "../include/utilities.h"

using namespace dealii;

/*Print 4th Order Symmetric Tensor*/
template <int dim>
void print_tensor(const SymmetricTensor<4,dim> &C){    
std::cout << "i" << "\t" <<"j" << "\t"<< "k" << "\t"<< "l" <<std::endl;
for(unsigned int i=0;i<dim;++i)
        for(unsigned int j=0;j<dim;++j)
                for(unsigned int k=0;k<dim;++k)
                        for(unsigned int l=0;l<dim;++l)
                        std::cout << i << "\t" << j << "\t"<< k << "\t"<< l << "\t\t" <<  C[i][j][k][l] << std::endl;       
        
}

/*Print 4th Order Tensor*/
template <int dim>
void print_tensor(const Tensor<4,dim> &C){
std::cout << "i" << "\t" <<"j" << "\t"<< "k" << "\t"<< "l" <<std::endl;
for(unsigned int i=0;i<dim;++i)
        for(unsigned int j=0;j<dim;++j)
                for(unsigned int k=0;k<dim;++k)
                        for(unsigned int l=0;l<dim;++l)
                        std::cout << i << "\t" << j << "\t"<< k << "\t"<< l << "\t\t" <<  C[i][j][k][l] << std::endl;       
        
}

/*Print 2nd Order Symmetric Tensor*/
template <int dim>
void print_tensor(const SymmetricTensor<2,dim> &C){
for(unsigned int i=0;i<dim;++i){
        for(unsigned int j=0;j<dim;++j)
               std::cout <<  C[i][j] <<"\t"; 
	std::cout<<std::endl;       
}
}

/*Print 2nd Order Tensor*/
template <int dim>
void print_tensor(const Tensor<2,dim> &C){
for(unsigned int i=0;i<dim;++i){
        for(unsigned int j=0;j<dim;++j)
                std::cout <<  C[i][j] <<"\t"; 
	std::cout<<std::endl;       
}
}
/*Returns a biaxial strain*/
template <int dim>
SymmetricTensor<2,dim> biaxial(){
        SymmetricTensor<2,dim> b;

        for(unsigned int i=0;i<dim;++i)
                for(unsigned int j=0;j<dim;++j){
                        if(i==j)
                        b[i][j]=0.2;
                        else
                        b[i][j]=0;
                }
return b;
}

/*Returns a uniaxial strain*/
template <int dim>
SymmetricTensor<2,dim> uniaxial(){
        SymmetricTensor<2,dim> u;
        for(unsigned int i=0;i<dim;++i)
                for(unsigned int j=0;j<dim;++j){
                        if(i==0 && j==0)
                        u[i][j]=0.2;
                        else
                        u[i][j]=0;
                }
return u;

}
/*To define the definition of sign funciton*/
double get_sign(double x)
{
        double sgn_x;

        if(std::fabs(x) < 1e-8)
		sgn_x = 0;
	else
		sgn_x = (x>0) ? 1:-1 ;
return sgn_x; 
}


template void print_tensor(const SymmetricTensor<4,2>&);
template void print_tensor(const SymmetricTensor<4,3>&);
template void print_tensor(const Tensor<4,2>&);
template void print_tensor(const Tensor<4,3>&);
template void print_tensor(const SymmetricTensor<2,2>&);
template void print_tensor(const SymmetricTensor<2,3>&);
template void print_tensor(const Tensor<2,2>&);
template void print_tensor(const Tensor<2,3>&);
template SymmetricTensor<2,2> biaxial();
template SymmetricTensor<2,3> biaxial();
template SymmetricTensor<2,2> uniaxial();
template SymmetricTensor<2,3> uniaxial();
