#include "../include/ElasticProblem.h"

using namespace dealii;

//Gives BigC without spectral decomposition as in LKM lecture notes
template <int dim>
SymmetricTensor<4, dim> get_const_BigC(const double lambda
                                                ,const double mu)
{

  SymmetricTensor<4, dim> C;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        for (unsigned int l = k; l < dim; ++l)
          C[i][j][k][l] = ((((i == k) && (j == l)) ? mu:0)
                            + (((i == l) && (j == k)) ? mu:0)
                             +((i==j) && (k==l)? lambda:0));

  return C;
}

//Print 4th Order Symmetric Tensor
template <int dim>
void print_C(SymmetricTensor<4,dim> &C){
for(unsigned int i=0;i<dim;++i)
	for(unsigned int j=0;j<dim;++j){
		for(unsigned int k=0;k<dim;++k)
			for(unsigned int l=0;l<dim;++l)
			std::cout<<C[i][j][k][l]<<"\t";
	std::cout<<std::endl;
	}
}

//Prints 4th Order Tensor
template <int dim>
void print_C(Tensor<4,dim> &C){
for(unsigned int i=0;i<dim;++i)
	for(unsigned int j=0;j<dim;++j){
		for(unsigned int k=0;k<dim;++k)
			for(unsigned int l=0;l<dim;++l)
			std::cout<<C[i][j][k][l]<<"\t";
	std::cout<<std::endl;
	}
}

//Prints strain 
template <int dim>
void print_eps(SymmetricTensor<2,dim> &eps){
for(unsigned int i=0;i<dim;++i){
	for(unsigned int j=0;j<dim;++j)
		std::cout<<eps[i][j]<<"\t";
	        std::cout<<std::endl;
	}
}

//Returns a biaxial strain
template <int dim>
SymmetricTensor<2,dim> biaxial(){
	SymmetricTensor<2,dim> b;
	
	for(unsigned int i=0;i<dim;++i)
		for(unsigned int j=0;j<dim;++j){
			if(i==j)
			b[i][j]=1;
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
                        u[i][j]=1;
                        else
                        u[i][j]=0;
                }
return u;

}


//Gives stress at particular quadrature point
template<int dim>
SymmetricTensor<2,dim> get_stress(const double lambda
				 ,const double mu
				 ,SymmetricTensor<2,dim> &eps){

	SymmetricTensor<2,dim> stress;
	double tr_eps = trace(eps);

	stress = lambda*tr_eps*unit_symmetric_tensor<dim>() + 2*mu*eps;

return stress; 
}

//Gives coefficient of first part of spectral BigC
double delsig_dellmbda(const double lambda
                      ,const double mu
		      ,unsigned int i
		      ,unsigned int j)
{
	double result = 0;

	if(i==j)
	result = lambda + 2*mu;
	else
	result = lambda;

return result;
}

//Gives part of coefficient after applying L'hospital rule (in second part of spectral BigC)
double delsig_dellmbda_b(const double lambda
                        ,const double mu)
{
        double result = 0;

        result =  2*mu + lambda;

return result;
}


//Gives part of +ve coefficient after applying L'hospital rule (in second part of +ve BigC)
template<int dim>
double delsig_dellmbda_b_plus(const double lambda
                             ,const double mu
                     	     ,SymmetricTensor<2,dim> & eps
			     ,double eps_eigenvalue)
{
        double result;
	double sgn_tr_eps, sgn_eps_eigenvalue;

	sgn_tr_eps = (trace(eps)>0) ? 1:0 ;
	sgn_eps_eigenvalue = (eps_eigenvalue>0) ? 1:0 ;

	result = 0.5*lambda*(1+sgn_tr_eps) + mu*(1+sgn_eps_eigenvalue);

return result;
}

//Gives part of -ve coefficient after applying L'hospital rule (in second part of -ve BigC)
template<int dim>
double delsig_dellmbda_b_minus(const double lambda
                             ,const double mu
                     	     ,SymmetricTensor<2,dim> & eps
			     ,double eps_eigenvalue)
{
        double result;
	double sgn_tr_eps, sgn_eps_eigenvalue;

	sgn_tr_eps = (trace(eps)>0) ? 1:0 ;
	sgn_eps_eigenvalue = (eps_eigenvalue>0) ? 1:0 ;

	result = 0.5*lambda*(1-sgn_tr_eps) + mu*(1-sgn_eps_eigenvalue);

return result;
}

/*Gives coefficient of first part of positive BigC*/
template<int dim>
double delsig_dellmbda_plus(const double lambda
                           ,const double mu
                     	   ,unsigned int i
                     	   ,unsigned int j
			   ,SymmetricTensor<2,dim> & eps
			   ,double eps_eigenvalue)
{
        double result;
	double sgn_tr_eps, sgn_eps_eigenvalue;

	sgn_tr_eps = (trace(eps)>0) ? 1:0 ;
	sgn_eps_eigenvalue = (eps_eigenvalue>0) ? 1:0 ;

	result = 0.5*lambda*(1+sgn_tr_eps) + ((i==j) ? mu*(1+sgn_eps_eigenvalue):0);

return result;
}

/*Gives coefficient of first part of negative BigC*/
template<int dim>
double delsig_dellmbda_minus(const double lambda
	                    ,const double mu
                     	    ,unsigned int i
                      	    ,unsigned int j
			    ,SymmetricTensor<2,dim> & eps
                            ,double eps_eigenvalue)
{
        double result;
        double sgn_tr_eps, sgn_eps_eigenvalue;

        sgn_tr_eps = (trace(eps)>0) ? 1:0 ;
        sgn_eps_eigenvalue = (eps_eigenvalue>0) ? 1:0 ;

	result = 0.5*lambda*(1-sgn_tr_eps) + ((i==j) ? mu*(1-sgn_eps_eigenvalue):0);

return result;
}

/*Gives stress eigenvalue to be used in coefficient of second part BigC*/
template<int dim>
double get_stress_eigenvalue(const double lambda
                      	    ,const double mu
                            ,SymmetricTensor<2,dim> & eps
			    ,double eps_eigenvalue)
{
	double result;
	double tr_eps = trace(eps);
	result = lambda*tr_eps + 2*mu*eps_eigenvalue;

return result;
}

/*Gives positive stress eigenvalue to be used in coefficient of second part positive BigC*/
template<int dim>
double get_stress_eigenvalue_plus(const double lambda
                                 ,const double mu
                      		 ,SymmetricTensor<2,dim> & eps
                        	 ,double eps_eigenvalue)
{
        double result;
        double tr_eps = trace(eps);

	double tr_eps_plus = (tr_eps>0) ? tr_eps:0 ;
	double eps_eigenvalue_plus = (eps_eigenvalue > 0) ? eps_eigenvalue:0 ;

        result = lambda*tr_eps_plus + 2*mu*eps_eigenvalue_plus;

return result;
}

/*Gives negative stress eigenvalue to be used in coefficient of second part negative BigC*/ 
template<int dim>
double get_stress_eigenvalue_minus(const double lambda
                      		  ,const double mu
                       		  ,SymmetricTensor<2,dim> & eps
                      		  ,double eps_eigenvalue)
{
        double result;
        double tr_eps = trace(eps);

        double tr_eps_minus = (tr_eps>0) ? 0:tr_eps;
        double eps_eigenvalue_minus = (eps_eigenvalue > 0) ? 0:eps_eigenvalue;

        result = lambda*tr_eps_minus + 2*mu*eps_eigenvalue_minus;

return result;
}

/*Gives BigC: both definitions works- with and without splitting into +ve & -ve*/
template <int dim>
SymmetricTensor<4,dim> get_BigC(const double lambda
                               ,const double mu
                               ,SymmetricTensor<2,dim> &eps)
{
	Tensor<4,dim> C_1;
	Tensor<4,dim> C_2;
	SymmetricTensor<4,dim> C_total;

	Tensor<4,dim> C_1_plus;
	Tensor<4,dim> C_1_minus;
	Tensor<4,dim> C_2_plus;
	Tensor<4,dim> C_2_minus;
	SymmetricTensor<4,dim> C_total_pm;

	Tensor<2,dim> Na_x_Na;
	Tensor<2,dim> Nb_x_Nb;
	Tensor<2,dim> Na_x_Nb;
	Tensor<2,dim> Nb_x_Na;

	std::array <std::pair< double, Tensor< 1, dim, double > >,std::integral_constant< int, dim >::value> eigen;
  	eigen = eigenvectors(eps,SymmetricTensorEigenvectorMethod::ql_implicit_shifts);


	//Calculating C_1:
	double scalar_1 = 0;
  	for (unsigned int a = 0;a < dim; ++a){
	 
		Na_x_Na = outer_product(eigen[a].second,eigen[a].second); 

		for (unsigned int b = 0;b < dim; ++b) {
	
			scalar_1 = delsig_dellmbda(lambda,mu,a,b);
			Nb_x_Nb = outer_product(eigen[b].second,eigen[b].second);
        
	  		C_1 +=  scalar_1*(outer_product(Na_x_Na,Nb_x_Nb));
      		}
	}


	//Calculating C_2:
	double scalar_2=0;
	double scalar_3=0;
        for (unsigned int a=0;a<dim;++a) 
                for (unsigned int b=0;b<dim;++b){
		scalar_2  = 0.5* ( get_stress_eigenvalue(lambda,mu,eps,eigen[b].first) - get_stress_eigenvalue(lambda,mu,eps,eigen[a].first))/(eigen[b].first - eigen[a].first);
		scalar_3 =  0.5* (delsig_dellmbda_b(lambda,mu) - delsig_dellmbda(lambda,mu,a,b));
		Na_x_Nb = outer_product(eigen[a].second,eigen[b].second);
		Nb_x_Na = outer_product(eigen[b].second,eigen[a].second);
			if(a!=b){
				if( !(std::fabs(eigen[a].first - eigen[b].first) < 1e-8)){
                        		C_2 += scalar_2 *( outer_product(Na_x_Nb,Na_x_Nb) + outer_product(Na_x_Nb,Nb_x_Na) );
				}
				else
					C_2 += scalar_3 *( outer_product(Na_x_Nb,Na_x_Nb) + outer_product(Na_x_Nb,Nb_x_Na) );
				}
       		}
 

	//Calculating total of C_1+C_2:
	for (unsigned int i = 0; i < dim; ++i)
		for (unsigned int j = i; j < dim; ++j)
			 for (unsigned int k = 0; k < dim; ++k)
				for (unsigned int l = k; l < dim; ++l)
		 		C_total[i][j][k][l] = C_1[i][j][k][l] + C_2[i][j][k][l];

	//Calculating C_1_plus:
	double scalar_4 = 0;
  	for (unsigned int a = 0;a < dim; ++a){
	 
		Na_x_Na = outer_product(eigen[a].second,eigen[a].second); 

		for (unsigned int b = 0;b < dim; ++b) {
	
			scalar_4 = delsig_dellmbda_plus(lambda,mu,a,b,eps,eigen[a].first);
			Nb_x_Nb = outer_product(eigen[b].second,eigen[b].second);
        
	  		C_1_plus +=  scalar_4*(outer_product(Na_x_Na,Nb_x_Nb));
      		}
	}

	//Calculating C_2_plus:
	double scalar_5=0;
	double scalar_6=0;
        for (unsigned int a=0;a<dim;++a) 
                for (unsigned int b=0;b<dim;++b){
		scalar_5  = 0.5* ( get_stress_eigenvalue_plus(lambda,mu,eps,eigen[b].first) - get_stress_eigenvalue_plus(lambda,mu,eps,eigen[a].first))/(eigen[b].first - eigen[a].first);
		scalar_6 =  0.5* (delsig_dellmbda_b_plus(lambda,mu,eps,eigen[b].first) - delsig_dellmbda_plus(lambda,mu,a,b,eps,eigen[a].first));
		Na_x_Nb = outer_product(eigen[a].second,eigen[b].second);
		Nb_x_Na = outer_product(eigen[b].second,eigen[a].second);
			if(a!=b){
				if( !(std::fabs(eigen[a].first - eigen[b].first) < 1e-8)){
                        		C_2_plus += scalar_5 *( outer_product(Na_x_Nb,Na_x_Nb) + outer_product(Na_x_Nb,Nb_x_Na) );
				}
				else
					C_2_plus += scalar_6 *( outer_product(Na_x_Nb,Na_x_Nb) + outer_product(Na_x_Nb,Nb_x_Na) );
				}
       		}
 



	//Calculating C_1_minus:
	double scalar_7 = 0;
  	for (unsigned int a = 0;a < dim; ++a){
	 
		Na_x_Na = outer_product(eigen[a].second,eigen[a].second); 

		for (unsigned int b = 0;b < dim; ++b) {
	
			scalar_7 = delsig_dellmbda_minus(lambda,mu,a,b,eps,eigen[a].first);
			Nb_x_Nb = outer_product(eigen[b].second,eigen[b].second);
        
	  		C_1_minus +=  scalar_7*(outer_product(Na_x_Na,Nb_x_Nb));
      		}
	}

	
	//Calculating C_2_mius:
	double scalar_8=0;
	double scalar_9=0;
        for (unsigned int a=0;a<dim;++a) 
                for (unsigned int b=0;b<dim;++b){
		scalar_8  = 0.5* ( get_stress_eigenvalue_minus(lambda,mu,eps,eigen[b].first) - get_stress_eigenvalue_minus(lambda,mu,eps,eigen[a].first))/(eigen[b].first - eigen[a].first);
		scalar_9 =  0.5* (delsig_dellmbda_b_minus(lambda,mu,eps,eigen[b].first) - delsig_dellmbda_minus(lambda,mu,a,b,eps,eigen[a].first));
		Na_x_Nb = outer_product(eigen[a].second,eigen[b].second);
		Nb_x_Na = outer_product(eigen[b].second,eigen[a].second);
			if(a!=b){
				if( !(std::fabs(eigen[a].first - eigen[b].first) < 1e-8)){
                        		C_2_minus += scalar_8 *( outer_product(Na_x_Nb,Na_x_Nb) + outer_product(Na_x_Nb,Nb_x_Na) );
				}
				else
					C_2_minus += scalar_9 *( outer_product(Na_x_Nb,Na_x_Nb) + outer_product(Na_x_Nb,Nb_x_Na) );
				}
       		}
 


	//Calculating total of C_1_plus + C_2_plus + C_1_minus + C_2_minus::
	for (unsigned int i = 0; i < dim; ++i)
		for (unsigned int j = i; j < dim; ++j)
			 for (unsigned int k = 0; k < dim; ++k)
				for (unsigned int l = k; l < dim; ++l)
		 		C_total_pm[i][j][k][l] = C_1_plus[i][j][k][l] + C_2_plus[i][j][k][l] + C_1_minus[i][j][k][l] + C_2_minus[i][j][k][l];


static bool once_2 = false;
if(!once_2){
std::cout<<"C_1:"<<std::endl;
print_C(C_1);
std::cout<<"C_2:"<<std::endl;
print_C(C_2);
once_2 =true;
}

return C_total;
//return C_total_pm;
}


template SymmetricTensor<4,2> get_const_BigC(double,double);
template SymmetricTensor<4,3> get_const_BigC(double,double);
template SymmetricTensor<4,2> get_BigC(double,double,SymmetricTensor<2,2>&);
template SymmetricTensor<4,3> get_BigC(double,double,SymmetricTensor<2,3>&);
template SymmetricTensor<2,2> get_stress(const double,const double,SymmetricTensor<2,2>&);
template SymmetricTensor<2,3> get_stress(const double,const double,SymmetricTensor<2,3>&);
template void print_C(SymmetricTensor<4,2>&);
template void print_C(SymmetricTensor<4,3>&);
template void print_C(Tensor<4,2>&);
template void print_C(Tensor<4,3>&);
template void print_eps(SymmetricTensor<2,2>&);
template void print_eps(SymmetricTensor<2,3>&);
template SymmetricTensor<2,2> biaxial();
template SymmetricTensor<2,3> biaxial();
template SymmetricTensor<2,2> uniaxial();
template SymmetricTensor<2,3> uniaxial();
