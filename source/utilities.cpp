#include "../include/utilities.h"

using namespace dealii;

//Print 4th Order Symmetric Tensor
template <int dim>
void print_tensor(SymmetricTensor<4,dim> &C){
for(unsigned int i=0;i<dim;++i)
        for(unsigned int j=0;j<dim;++j)
                for(unsigned int k=0;k<dim;++k)
                        for(unsigned int l=0;l<dim;++l)
                        std::cout << i << "\t" << j << "\t"<< k << "\t"<< l << "\t\t" <<  C[i][j][k][l] << std::endl;       
        
}

//Print 4th Order Tensor
template <int dim>
void print_tensor(Tensor<4,dim> &C){
for(unsigned int i=0;i<dim;++i)
        for(unsigned int j=0;j<dim;++j)
                for(unsigned int k=0;k<dim;++k)
                        for(unsigned int l=0;l<dim;++l)
                        std::cout << i << "\t" << j << "\t"<< k << "\t"<< l << "\t\t" <<  C[i][j][k][l] << std::endl;       
        
}

//Print 2nd Order Symmetric Tensor
template <int dim>
void print_tensor(SymmetricTensor<2,dim> &C){
for(unsigned int i=0;i<dim;++i)
        for(unsigned int j=0;j<dim;++j){
               std::cout << i << "\t" << j << "\t\t" <<  C[i][j] << std::endl;       
        }
}

//Print 2nd Order Tensor
template <int dim>
void print_tensor(Tensor<2,dim> &C){
for(unsigned int i=0;i<dim;++i)
        for(unsigned int j=0;j<dim;++j){
               std::cout << i << "\t" << j << "\t\t" <<  C[i][j] << std::endl;       
        }
}



template void print_tensor(SymmetricTensor<4,2>&);
template void print_tensor(SymmetricTensor<4,3>&);
template void print_tensor(Tensor<4,2>&);
template void print_tensor(Tensor<4,3>&);
template void print_tensor(SymmetricTensor<2,2>&);
template void print_tensor(SymmetricTensor<2,3>&);
template void print_tensor(Tensor<2,2>&);
template void print_tensor(Tensor<2,3>&);

