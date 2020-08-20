#include "../include/utilities.h"

using namespace dealii;

//Print 4th Order Symmetric Tensor
template <int dim>
void print_tensor(const SymmetricTensor<4,dim> &C){    
std::cout << "i" << "\t" <<"j" << "\t"<< "k" << "\t"<< "l" <<std::endl;
for(unsigned int i=0;i<dim;++i)
        for(unsigned int j=0;j<dim;++j)
                for(unsigned int k=0;k<dim;++k)
                        for(unsigned int l=0;l<dim;++l)
                        std::cout << i << "\t" << j << "\t"<< k << "\t"<< l << "\t\t" <<  C[i][j][k][l] << std::endl;       
        
}

//Print 4th Order Tensor
template <int dim>
void print_tensor(const Tensor<4,dim> &C){
std::cout << "i" << "\t" <<"j" << "\t"<< "k" << "\t"<< "l" <<std::endl;
for(unsigned int i=0;i<dim;++i)
        for(unsigned int j=0;j<dim;++j)
                for(unsigned int k=0;k<dim;++k)
                        for(unsigned int l=0;l<dim;++l)
                        std::cout << i << "\t" << j << "\t"<< k << "\t"<< l << "\t\t" <<  C[i][j][k][l] << std::endl;       
        
}

//Print 2nd Order Symmetric Tensor
template <int dim>
void print_tensor(const SymmetricTensor<2,dim> &C){
for(unsigned int i=0;i<dim;++i){
        for(unsigned int j=0;j<dim;++j)
               std::cout <<  C[i][j] <<"\t"; 
	std::cout<<std::endl;       
}
}

//Print 2nd Order Tensor
template <int dim>
void print_tensor(const Tensor<2,dim> &C){
for(unsigned int i=0;i<dim;++i){
        for(unsigned int j=0;j<dim;++j)
                std::cout <<  C[i][j] <<"\t"; 
	std::cout<<std::endl;       
}
}

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

//Returns a uniaxial strain
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

double get_sign(double x)
{
        double sgn_x;

        if(std::fabs(x) < 1e-8)
		sgn_x = 0;
	else
		sgn_x = (x>0) ? 1:-1 ;
return sgn_x; 
}

template<int dim>
void comparison(const double lambda,const double mu,SymmetricTensor<2,dim> &dummy)
{
  SymmetricTensor<2,dim> eps0,eps1,eps2,eps3,eps4,eps5;
   
  for(unsigned int i=0;i<dim;++i)
	  for(unsigned int j=0;j<dim;++j)
	  {
		eps0[i][j]=0;
	  }	
	
  for(unsigned int i=0;i<dim;++i)

	  for(unsigned int j=0;j<dim;++j)
	  {
		  if(i==0 && j==0)
			  eps1[i][j]=-0.1;
		  else if(i==(dim-1) && j==(dim-1))
			  eps1[i][j]=0.2;
		  else
			  eps1[i][j]=0;
	  }

 for(unsigned int i=0;i<dim;++i)
	  for(unsigned int j=0;j<dim;++j)
	  {
		  if(i != j)
			  eps2[i][j]=-0.1;
		  else
			  eps2[i][j]=0;
	  }

 
  for(unsigned int i=0;i<dim;++i)
	  for(unsigned int j=0;j<dim;++j)
	  {
		  if(i==0 && j==0)
			  eps3[i][j]=0.2;
		  else if(i==(dim-1) && j==(dim-1))
			  eps3[i][j]=0.1;
		  else
			  eps3[i][j]=0.05;
	  }

 for(unsigned int i=0;i<dim;++i)
	  for(unsigned int j=0;j<dim;++j)
	  {
		  if(i != j)
			  eps4[i][j]=0.1;
		  else
			  eps4[i][j]=0;
	  }


 for(unsigned int i=0;i<dim;++i)
                for(unsigned int j=0;j<dim;++j){
                        if(i==j)
                        eps5[i][j]=-0.1;
                        else
                        eps5[i][j]=0;
                }
//Printing

std::cout<<std::endl;
std::cout<<"Epsilon:"<<std::endl;
print_tensor(eps0);
std::cout<<"BigC_plus:"<<std::endl;
print_tensor( get_BigC_plus(lambda,mu,eps0) );
std::cout<<"Norm:"<<get_BigC_plus(lambda,mu,eps0).norm()<<std::endl;
std::cout<<std::endl;

std::cout<<"Epsilon:"<<std::endl;
print_tensor(eps1);
std::cout<<"BigC_plus:"<<std::endl;
print_tensor( get_BigC_plus(lambda,mu,eps1) );
std::cout<<"Norm:"<<get_BigC_plus(lambda,mu,eps1).norm()<<std::endl;
std::cout<<std::endl;


std::cout<<std::endl;
std::cout<<"Epsilon:"<<std::endl;
print_tensor(eps2);
std::cout<<"BigC_plus:"<<std::endl;
print_tensor( get_BigC_plus(lambda,mu,eps2) );
std::cout<<"Norm:"<<get_BigC_plus(lambda,mu,eps2).norm()<<std::endl;
std::cout<<std::endl;


std::cout<<std::endl;
std::cout<<"Epsilon:"<<std::endl;
print_tensor(eps3);
std::cout<<"BigC_plus:"<<std::endl;
print_tensor( get_BigC_plus(lambda,mu,eps3) );
std::cout<<"Norm:"<<get_BigC_plus(lambda,mu,eps3).norm()<<std::endl;
std::cout<<std::endl;


std::cout<<std::endl;
std::cout<<"Epsilon:"<<std::endl;
print_tensor(eps4);
std::cout<<"BigC_plus:"<<std::endl;
print_tensor( get_BigC_plus(lambda,mu,eps4) );
std::cout<<"Norm:"<<get_BigC_plus(lambda,mu,eps4).norm()<<std::endl;
std::cout<<std::endl;


std::cout<<"Epsilon:"<<std::endl;
print_tensor(eps5);
std::cout<<"BigC_plus:"<<std::endl;
print_tensor( get_BigC_plus(lambda,mu,eps5) );
std::cout<<"Norm:"<<get_BigC_plus(lambda,mu,eps5).norm()<<std::endl;
std::cout<<std::endl;


std::cout<<"------------------------------------------------"<<std::endl;
std::cout<<std::endl;
std::cout<<"Epsilon:"<<std::endl;
print_tensor(eps0);
std::cout<<"BigC_minus:"<<std::endl;
print_tensor( get_BigC_minus(lambda,mu,eps0) );
std::cout<<"Norm:"<<get_BigC_minus(lambda,mu,eps0).norm()<<std::endl;
std::cout<<std::endl;


std::cout<<"Epsilon:"<<std::endl;
print_tensor(eps1);
std::cout<<"BigC_minus:"<<std::endl;
print_tensor( get_BigC_minus(lambda,mu,eps1) );
std::cout<<"Norm:"<<get_BigC_minus(lambda,mu,eps1).norm()<<std::endl;
std::cout<<std::endl;


std::cout<<std::endl;
std::cout<<"Epsilon:"<<std::endl;
print_tensor(eps2);
std::cout<<"BigC_minus:"<<std::endl;
print_tensor( get_BigC_minus(lambda,mu,eps2) );
std::cout<<"Norm:"<<get_BigC_minus(lambda,mu,eps2).norm()<<std::endl;
std::cout<<std::endl;


std::cout<<std::endl;
std::cout<<"Epsilon:"<<std::endl;
print_tensor(eps3);
std::cout<<"BigC_minus:"<<std::endl;
print_tensor( get_BigC_minus(lambda,mu,eps3) );
std::cout<<"Norm:"<<get_BigC_minus(lambda,mu,eps3).norm()<<std::endl;
std::cout<<std::endl;


std::cout<<std::endl;
std::cout<<"Epsilon:"<<std::endl;
print_tensor(eps4);
std::cout<<"BigC_minus:"<<std::endl;
print_tensor( get_BigC_minus(lambda,mu,eps4) );
std::cout<<"Norm:"<<get_BigC_minus(lambda,mu,eps4).norm()<<std::endl;
std::cout<<std::endl;


std::cout<<"Epsilon:"<<std::endl;
print_tensor(eps5);
std::cout<<"BigC_minus:"<<std::endl;
print_tensor( get_BigC_minus(lambda,mu,eps5) );
std::cout<<"Norm:"<<get_BigC_minus(lambda,mu,eps5).norm()<<std::endl;
std::cout<<std::endl;



std::cout<<"------------------------------------------------"<<std::endl;
std::cout<<"Const_BigC"<<std::endl;
print_tensor(get_const_BigC(lambda,mu,dummy));
std::cout<<"Norm:"<<get_const_BigC(lambda,mu,dummy).norm()<<std::endl;
std::cout<<std::endl;


std::cout<<std::endl;
std::cout<<"Epsilon:"<<std::endl;
print_tensor(eps0);
std::cout<<"BigC_total_pm:"<<std::endl;
print_tensor( get_BigC(lambda,mu,eps0) );
std::cout<<"Norm:"<<get_BigC(lambda,mu,eps0).norm()<<std::endl;
std::cout<<std::endl;


std::cout<<"Epsilon:"<<std::endl;
print_tensor(eps1);
std::cout<<"BigC_total_pm:"<<std::endl;
print_tensor( get_BigC(lambda,mu,eps1) );
std::cout<<"Norm:"<<get_BigC(lambda,mu,eps1).norm()<<std::endl;
std::cout<<std::endl;


std::cout<<std::endl;
std::cout<<"Epsilon:"<<std::endl;
print_tensor(eps2);
std::cout<<"BigC_total_pm:"<<std::endl;
print_tensor( get_BigC(lambda,mu,eps2) );
std::cout<<"Norm:"<<get_BigC(lambda,mu,eps2).norm()<<std::endl;
std::cout<<std::endl;


std::cout<<std::endl;
std::cout<<"Epsilon:"<<std::endl;
print_tensor(eps3);
std::cout<<"BigC_total_pm:"<<std::endl;
print_tensor( get_BigC(lambda,mu,eps3) );
std::cout<<"Norm:"<<get_BigC(lambda,mu,eps3).norm()<<std::endl;
std::cout<<std::endl;


std::cout<<std::endl;
std::cout<<"Epsilon:"<<std::endl;
print_tensor(eps4);
std::cout<<"BigC_total_pm:"<<std::endl;
print_tensor( get_BigC(lambda,mu,eps4) );
std::cout<<"Norm:"<<get_BigC(lambda,mu,eps4).norm()<<std::endl;
std::cout<<std::endl;


std::cout<<"Epsilon:"<<std::endl;
print_tensor(eps5);
std::cout<<"BigC_total_pm:"<<std::endl;
print_tensor( get_BigC(lambda,mu,eps5) );
std::cout<<"Norm:"<<get_BigC(lambda,mu,eps5).norm()<<std::endl;
std::cout<<std::endl;

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
template void comparison(const double ,const double ,SymmetricTensor<2,2> &);
template void comparison(const double ,const double ,SymmetricTensor<2,3> &);
