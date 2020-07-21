#include "../include/utilities.h"

using namespace dealii;

//Gives BigC without spectral decomposition as in LKM lecture notes
template <int dim>
SymmetricTensor<4, dim> get_const_BigC(const double lambda
                                      ,const double mu
				      ,SymmetricTensor<2,dim> &dummy)
{

  SymmetricTensor<4, dim> C;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        for (unsigned int l = k; l < dim; ++l)
          C[i][j][k][l] = ((((i == k) && (j == l)) ? mu:0)
                            + (((i == l) && (j == k)) ? mu:0)
                             +((i==j) && (k==l)? lambda:0));
/*Printing the tensor
  static bool once_1=false;
  if(!once_1){
  std::cout<<"Const Big_C:"<<std::endl;
  print_tensor(C);
  once_1 =true;
  }
*/ 
  return C;
}

//Gives stress_plus at particular quadrature point
template<int dim>
SymmetricTensor<2,dim> get_stress_plus(const double lambda
				      ,const double mu
				      ,SymmetricTensor<2,dim> &eps){

	SymmetricTensor<2,dim> stress_plus,eps_plus;
	double tr_eps = trace(eps);
	double tr_eps_plus = (tr_eps>0) ? tr_eps:0;

	Tensor<2,dim> A;
	double scalar,scalar_plus;
	Tensor<2,dim> Na_x_Na;
	//Gives an array of dim pair of eigenvalues and eigenvectors of eps 
	std::array <std::pair< double, Tensor< 1, dim, double > >,std::integral_constant< int, dim >::value> eigen;
  	eigen = eigenvectors(eps,SymmetricTensorEigenvectorMethod::ql_implicit_shifts);

  	for (unsigned int a = 0;a < dim; ++a){
	 	Na_x_Na = outer_product(eigen[a].second,eigen[a].second); 
		scalar = eigen[a].first;
		scalar_plus = (scalar>0) ? scalar:0;
		A += scalar_plus*Na_x_Na;
	}

	for (unsigned int i = 0; i < dim; ++i)
		for (unsigned int j = 0; j < dim; ++j)
	 		eps_plus[i][j] = A[i][j];


	stress_plus = lambda*tr_eps_plus*unit_symmetric_tensor<dim>() + 2*mu*eps_plus;
/*
//Printing///////////////////  
  static bool once_1=false;
  if(!once_1){
  std::cout<<"Stress_plus:"<<std::endl;
  print_tensor(stress_plus);
  once_1 =true;
  }
//////////////////////////
*/
return stress_plus; 
}

//Gives stress_minus at particular quadrature point
template<int dim>
SymmetricTensor<2,dim> get_stress_minus(const double lambda
				      ,const double mu
				      ,SymmetricTensor<2,dim> &eps){

	SymmetricTensor<2,dim> stress_minus,eps_minus;
	double tr_eps = trace(eps);
	double tr_eps_minus = (tr_eps>0) ? 0:tr_eps;

	Tensor<2,dim> A;
	double scalar,scalar_minus;
	Tensor<2,dim> Na_x_Na;
	//Gives an array of dim pair of eigenvalues and eigenvectors of eps 
	std::array <std::pair< double, Tensor< 1, dim, double > >,std::integral_constant< int, dim >::value> eigen;
  	eigen = eigenvectors(eps,SymmetricTensorEigenvectorMethod::ql_implicit_shifts);

  	for (unsigned int a = 0;a < dim; ++a){
	 	Na_x_Na = outer_product(eigen[a].second,eigen[a].second); 
		scalar = eigen[a].first;
		scalar_minus = (scalar>0) ? 0:scalar;
		A += scalar_minus*Na_x_Na;
	}

	for (unsigned int i = 0; i < dim; ++i)
		for (unsigned int j = 0; j < dim; ++j)
	 		eps_minus[i][j] = A[i][j];


	stress_minus = lambda*tr_eps_minus*unit_symmetric_tensor<dim>() + 2*mu*eps_minus;
/*
//Printing///////////////////  
  static bool once_1=false;
  if(!once_1){
  std::cout<<"Stress_minus:"<<std::endl;
  print_tensor(stress_minus);
  once_1 =true;
  }
//////////////////////////
*/
return stress_minus; 
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
/*
//Gives an array of stress eigenvalues
template <int dim>
std::array<double,std::integral_constant< int, dim >::value> get_stress_eigenval(const double lambda
				      						,const double mu
				      						,SymmetricTensor<2,dim> &eps
										,std::array<double,std::integral_constant< int, dim >::value> &eps_eigenval){

	std::array<double,std::integral_constant< int, dim >::value> stress_eigenval;
	double tr_eps = trace(eps);

	for(unsigned int i=0;i<dim;++i)
		stress_eigenval[i] = lambda*tr_eps + 2*mu*eps_eigenval[i];
return stress_eigenval;
}
*/
//Gives an array of +ve stress eigenvalues
template <int dim>
std::array<double,std::integral_constant< int, dim >::value> get_stress_eigenval_plus(const double lambda
				      						,const double mu
				      						,SymmetricTensor<2,dim> &eps
										,std::array<double,std::integral_constant< int, dim >::value> &eps_eigenval){

	std::array<double,std::integral_constant< int, dim >::value> stress_eigenval_plus;
	double tr_eps = trace(eps);
	double tr_eps_plus = (tr_eps>0) ? tr_eps:0 ;

	std::array<double,std::integral_constant< int, dim >::value> eps_eigenval_plus;

	for(unsigned int i=0;i<dim;++i)
		eps_eigenval_plus[i] = (eps_eigenval[i] > 0) ? eps_eigenval[i]:0 ;

	for(unsigned int i=0;i<dim;++i)
		stress_eigenval_plus[i] = lambda*tr_eps_plus + 2*mu*eps_eigenval_plus[i];

return stress_eigenval_plus;
}

//Gives an array of -ve stress eigenvalues
template <int dim>
std::array<double,std::integral_constant< int, dim >::value> get_stress_eigenval_minus(const double lambda
				      						,const double mu
				      						,SymmetricTensor<2,dim> &eps
										,std::array<double,std::integral_constant< int, dim >::value> &eps_eigenval){

	std::array<double,std::integral_constant< int, dim >::value> stress_eigenval_minus;
	double tr_eps = trace(eps);
	double tr_eps_minus = (tr_eps>0) ? 0:tr_eps ;

	std::array<double,std::integral_constant< int, dim >::value> eps_eigenval_minus;

	for(unsigned int i=0;i<dim;++i)
		eps_eigenval_minus[i] = (eps_eigenval[i] > 0) ? 0:eps_eigenval[i] ;

	for(unsigned int i=0;i<dim;++i)
		stress_eigenval_minus[i] = lambda*tr_eps_minus + 2*mu*eps_eigenval_minus[i];

return stress_eigenval_minus;
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

template <int dim>
SymmetricTensor<4,dim> get_BigC_plus(const double lambda
                                            ,const double mu
                                            ,dealii::SymmetricTensor<2,dim>& eps)
{

	Tensor<4,dim> C_1_plus;
	Tensor<4,dim> C_2_plus;
	SymmetricTensor<4,dim> C_total_plus;

	Tensor<2,dim> Na_x_Na;
	Tensor<2,dim> Nb_x_Nb;
	Tensor<2,dim> Na_x_Nb;
	Tensor<2,dim> Nb_x_Na;


	//Gives an array of dim pair of eigenvalues and eigenvectors of eps 
	std::array <std::pair< double, Tensor< 1, dim, double > >,std::integral_constant< int, dim >::value> eigen;
  	eigen = eigenvectors(eps,SymmetricTensorEigenvectorMethod::ql_implicit_shifts);


	//Gives an array of eps eigenvalues
	std::array <double,std::integral_constant< int, dim >::value> eps_eigenval;
	eps_eigenval = eigenvalues(eps);

	//Calculates an array of +ve stress eigenvalues
	std::array <double,std::integral_constant< int, dim >::value> stress_eigenval_plus;
	stress_eigenval_plus = get_stress_eigenval_plus(lambda,mu,eps,eps_eigenval); 
	

	//Calculating C_1_plus:
	double scalar_10 = 0;
  	for (unsigned int a = 0;a < dim; ++a){
	 
		Na_x_Na = outer_product(eigen[a].second,eigen[a].second); 

		for (unsigned int b = 0;b < dim; ++b) {
	
			scalar_10 = delsig_dellmbda_plus(lambda,mu,a,b,eps,eigen[a].first);
			Nb_x_Nb = outer_product(eigen[b].second,eigen[b].second);
        
	  		C_1_plus +=  scalar_10*(outer_product(Na_x_Na,Nb_x_Nb));
      		}
	}

	//Calculating C_2_plus:
	double scalar_11=0;
	double scalar_12=0;
        for (unsigned int a=0;a<dim;++a) 
                for (unsigned int b=0;b<dim;++b){
		scalar_12 = 0.5* (delsig_dellmbda_b_plus(lambda,mu,eps,eigen[b].first) - delsig_dellmbda_plus(lambda,mu,a,b,eps,eigen[a].first));
		Na_x_Nb = outer_product(eigen[a].second,eigen[b].second);
		Nb_x_Na = outer_product(eigen[b].second,eigen[a].second);
			if(a!=b){
		
				scalar_11 = 0.5* ( stress_eigenval_plus[b] - stress_eigenval_plus[a] )/(eigen[b].first - eigen[a].first);

				if( !(std::fabs(eigen[a].first - eigen[b].first) < 1e-8)){
                        		C_2_plus += scalar_11 *( outer_product(Na_x_Nb,Na_x_Nb) + outer_product(Na_x_Nb,Nb_x_Na) );
				}
				else
					C_2_plus += scalar_12 *( outer_product(Na_x_Nb,Na_x_Nb) + outer_product(Na_x_Nb,Nb_x_Na) );
				}
       		}
 


	//Calculating total of C_1_plus + C_2_plus :
	for (unsigned int i = 0; i < dim; ++i)
		for (unsigned int j = i; j < dim; ++j)
			 for (unsigned int k = 0; k < dim; ++k)
				for (unsigned int l = k; l < dim; ++l)
		 		C_total_plus[i][j][k][l] = C_1_plus[i][j][k][l] + C_2_plus[i][j][k][l];
/*
///Printing//////////////////
  static bool once_1=false;
  if(!once_1){
  std::cout<<"Big_C_plus:"<<std::endl;
  print_tensor(C_total_plus);
  once_1 =true;
  }
//////////////////////////
*/
return C_total_plus;
}


template <int dim>
SymmetricTensor<4,dim> get_BigC_minus(const double lambda
                                             ,const double mu
                                             ,dealii::SymmetricTensor<2,dim>& eps)
{
	Tensor<4,dim> C_1_minus;
	Tensor<4,dim> C_2_minus;
	SymmetricTensor<4,dim> C_total_minus;

	Tensor<2,dim> Na_x_Na;
	Tensor<2,dim> Nb_x_Nb;
	Tensor<2,dim> Na_x_Nb;
	Tensor<2,dim> Nb_x_Na;

	//Gives an array of dim pair of eigenvalues and eigenvectors of eps 
	std::array <std::pair< double, Tensor< 1, dim, double > >,std::integral_constant< int, dim >::value> eigen;
  	eigen = eigenvectors(eps,SymmetricTensorEigenvectorMethod::ql_implicit_shifts);

	//Gives an array of eps eigenvalues
	std::array <double,std::integral_constant< int, dim >::value> eps_eigenval;
	eps_eigenval = eigenvalues(eps);


	//Calculates an array of -ve stress eigenvalues
	std::array <double,std::integral_constant< int, dim >::value> stress_eigenval_minus;
	stress_eigenval_minus = get_stress_eigenval_minus(lambda,mu,eps,eps_eigenval); 
	

	//Calculating C_1_minus:
	double scalar_13 = 0;
  	for (unsigned int a = 0;a < dim; ++a){
	 
		Na_x_Na = outer_product(eigen[a].second,eigen[a].second); 

		for (unsigned int b = 0;b < dim; ++b) {
	
			scalar_13 = delsig_dellmbda_minus(lambda,mu,a,b,eps,eigen[a].first);
			Nb_x_Nb = outer_product(eigen[b].second,eigen[b].second);
        
	  		C_1_minus +=  scalar_13*(outer_product(Na_x_Na,Nb_x_Nb));
      		}
	}

	
	//Calculating C_2_mius:
	double scalar_14=0;
	double scalar_15=0;
        for (unsigned int a=0;a<dim;++a) 
                for (unsigned int b=0;b<dim;++b){
		scalar_15 =  0.5* (delsig_dellmbda_b_minus(lambda,mu,eps,eigen[b].first) - delsig_dellmbda_minus(lambda,mu,a,b,eps,eigen[a].first));
		Na_x_Nb = outer_product(eigen[a].second,eigen[b].second);
		Nb_x_Na = outer_product(eigen[b].second,eigen[a].second);
			if(a!=b){
		
				scalar_14  = 0.5* ( stress_eigenval_minus[b] - stress_eigenval_minus[a] )/(eigen[b].first - eigen[a].first);

				if( !(std::fabs(eigen[a].first - eigen[b].first) < 1e-8)){
                        		C_2_minus += scalar_14 *( outer_product(Na_x_Nb,Na_x_Nb) + outer_product(Na_x_Nb,Nb_x_Na) );
				}
				else
					C_2_minus += scalar_15 *( outer_product(Na_x_Nb,Na_x_Nb) + outer_product(Na_x_Nb,Nb_x_Na) );
				}
       		}
 


	//Calculating total of C_1_minus + C_2_minus::
	for (unsigned int i = 0; i < dim; ++i)
		for (unsigned int j = i; j < dim; ++j)
			 for (unsigned int k = 0; k < dim; ++k)
				for (unsigned int l = k; l < dim; ++l)
		 		C_total_minus[i][j][k][l] = C_1_minus[i][j][k][l] + C_2_minus[i][j][k][l];
/*
/////Printing/////////////
  static bool once_1=false;
  if(!once_1){
  std::cout<<"Big_C_minus:"<<std::endl;
  print_tensor(C_total_minus);
  once_1 =true;
  }
/////////////////////////
*/
return C_total_minus;
}



/*Gives BigC: both definitions works- with and without splitting into +ve & -ve*/
template <int dim>
SymmetricTensor<4,dim> get_BigC(const double lambda
                               ,const double mu
                               ,SymmetricTensor<2,dim> &eps)
{
/*	Tensor<4,dim> C_1;
	Tensor<4,dim> C_2;
	SymmetricTensor<4,dim> C_total;
*/
	Tensor<4,dim> C_1_plus;
	Tensor<4,dim> C_1_minus;
	Tensor<4,dim> C_2_plus;
	Tensor<4,dim> C_2_minus;
	SymmetricTensor<4,dim> C_total_pm;

	Tensor<2,dim> Na_x_Na;
	Tensor<2,dim> Nb_x_Nb;
	Tensor<2,dim> Na_x_Nb;
	Tensor<2,dim> Nb_x_Na;

	//Gives an array of dim pair of eigenvalues and eigenvectors of eps 
	std::array <std::pair< double, Tensor< 1, dim, double > >,std::integral_constant< int, dim >::value> eigen;
  	eigen = eigenvectors(eps,SymmetricTensorEigenvectorMethod::ql_implicit_shifts);
/*
	//Gives an array of eps eigenvalues
	std::array <double,std::integral_constant< int, dim >::value> eps_eigenval;
	eps_eigenval = eigenvalues(eps);

	//Calculates an array of stress eigenvalues
	std::array <double,std::integral_constant< int, dim >::value> stress_eigenval;
	stress_eigenval = get_stress_eigenval(lambda,mu,eps,eps_eigenval); 
	
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
		scalar_2  = 0.5* ( stress_eigenval[b] - stress_eigenval[a] )/(eigen[b].first - eigen[a].first);
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


	//Calculates an array of +ve stress eigenvalues
	std::array <double,std::integral_constant< int, dim >::value> stress_eigenval_plus;
	stress_eigenval_plus = get_stress_eigenval_plus(lambda,mu,eps,eps_eigenval); 
*/	

	//Testing PhaseField get_stress_plus function

	std::array <std::pair< double, Tensor< 1, dim, double > >,std::integral_constant< int, dim >::value> eigen_stress_p;
  	eigen_stress_p = eigenvectors(get_stress_plus(lambda,mu,eps),SymmetricTensorEigenvectorMethod::ql_implicit_shifts);

	std::array <double,std::integral_constant< int, dim >::value> stress_eigenval_plus;
	for(unsigned int c=0;c<stress_eigenval_plus.size();++c)
	stress_eigenval_plus[c] = eigen_stress_p[c].first; 
	

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
//		scalar_5  = 0.5* ( stress_eigenval_plus[b] - stress_eigenval_plus[a] )/(eigen[b].first - eigen[a].first);//(a==b)value nvrused
//		std::cout<<"scalar_5"<<a<<b<<" "<<scalar_5<<std::endl;
//		std::cout<<"eigen[b].first, eigen[a].first:"<<eigen[b].first<<" "<<eigen[a].first<<std::endl;
		scalar_6 =  0.5* (delsig_dellmbda_b_plus(lambda,mu,eps,eigen[b].first) - delsig_dellmbda_plus(lambda,mu,a,b,eps,eigen[a].first));
		Na_x_Nb = outer_product(eigen[a].second,eigen[b].second);
		Nb_x_Na = outer_product(eigen[b].second,eigen[a].second);
			if(a!=b){
			
				scalar_5  = 0.5* ( stress_eigenval_plus[b] - stress_eigenval_plus[a] )/(eigen[b].first - eigen[a].first);
			
				if( !(std::fabs(eigen[a].first - eigen[b].first) < 1e-8)){
                        		C_2_plus += scalar_5 *( outer_product(Na_x_Nb,Na_x_Nb) + outer_product(Na_x_Nb,Nb_x_Na) );
				}
				else
					C_2_plus += scalar_6 *( outer_product(Na_x_Nb,Na_x_Nb) + outer_product(Na_x_Nb,Nb_x_Na) );
				}
       		}
 

	//Calculates an array of -ve stress eigenvalues
//	std::array <double,std::integral_constant< int, dim >::value> stress_eigenval_minus;
//	stress_eigenval_minus = get_stress_eigenval_minus(lambda,mu,eps,eps_eigenval); 
	
	//Testing PhaseField get_stress_minus function

	std::array <std::pair< double, Tensor< 1, dim, double > >,std::integral_constant< int, dim >::value> eigen_stress_m;
  	eigen_stress_m = eigenvectors(get_stress_minus(lambda,mu,eps),SymmetricTensorEigenvectorMethod::ql_implicit_shifts);

	std::array <double,std::integral_constant< int, dim >::value> stress_eigenval_minus;
	for(unsigned int c=0;c<stress_eigenval_minus.size();++c)
	stress_eigenval_minus[c] = eigen_stress_m[c].first; 
	

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

		scalar_9 =  0.5* (delsig_dellmbda_b_minus(lambda,mu,eps,eigen[b].first) - delsig_dellmbda_minus(lambda,mu,a,b,eps,eigen[a].first));
		Na_x_Nb = outer_product(eigen[a].second,eigen[b].second);
		Nb_x_Na = outer_product(eigen[b].second,eigen[a].second);
			if(a!=b){

				scalar_8  = 0.5* ( stress_eigenval_minus[b] - stress_eigenval_minus[a] )/(eigen[b].first - eigen[a].first);

				if( !(std::fabs(eigen[a].first - eigen[b].first) < 1e-8)){
                        		C_2_minus += scalar_8 *( outer_product(Na_x_Nb,Na_x_Nb) + outer_product(Na_x_Nb,Nb_x_Na) );
				}
				else
					C_2_minus += scalar_9 *( outer_product(Na_x_Nb,Na_x_Nb) + outer_product(Na_x_Nb,Nb_x_Na) );
				}
       		}
 
	//
	unsigned int g=2;

	//Calculating total of C_1_plus + C_2_plus + C_1_minus + C_2_minus::
	for (unsigned int i = 0; i < dim; ++i)
		for (unsigned int j = i; j < dim; ++j)
			 for (unsigned int k = 0; k < dim; ++k)
				for (unsigned int l = k; l < dim; ++l)
		 		C_total_pm[i][j][k][l] = g*( C_1_plus[i][j][k][l] + C_2_plus[i][j][k][l]) + C_1_minus[i][j][k][l] + C_2_minus[i][j][k][l];

/*Printing Different 4th order tensors:
static bool once_2 = false;
if(!once_2){
std::cout<<"C_1:"<<std::endl;
print_tensor(C_1);
std::cout<<"C_2:"<<std::endl;
print_tensor(C_2);
std::cout<<"eps bigC (C_1+C_2):"<<std::endl;
print_tensor(C_total);
once_2 =true;
}
*/

//return C_total;
return C_total_pm;

//Testing Phasefield get_BigC_plus/minus
//SymmetricTensor<4,dim> tmp1 = get_BigC_minus(lambda,mu,eps);
//SymmetricTensor<4,dim> tmp2 = get_BigC_plus(lambda,mu,eps);
//return tmp1+tmp2;
}

template SymmetricTensor<4,2> get_const_BigC(double,double,SymmetricTensor<2,2>&);
template SymmetricTensor<4,3> get_const_BigC(double,double,SymmetricTensor<2,3>&);
template SymmetricTensor<4,2> get_BigC(double,double,SymmetricTensor<2,2>&);
template SymmetricTensor<4,3> get_BigC(double,double,SymmetricTensor<2,3>&);
template SymmetricTensor<2,2> get_stress(const double,const double,SymmetricTensor<2,2>&);
template SymmetricTensor<2,3> get_stress(const double,const double,SymmetricTensor<2,3>&);
template SymmetricTensor<4,2> get_BigC_plus(double,double,SymmetricTensor<2,2>&);
template SymmetricTensor<4,3> get_BigC_plus(double,double,SymmetricTensor<2,3>&);
template SymmetricTensor<4,2> get_BigC_minus(double,double,SymmetricTensor<2,2>&);
template SymmetricTensor<4,3> get_BigC_minus(double,double,SymmetricTensor<2,3>&);
template SymmetricTensor<2,2> get_stress_plus(const double,const double,SymmetricTensor<2,2>&);
template SymmetricTensor<2,3> get_stress_plus(const double,const double,SymmetricTensor<2,3>&);
template SymmetricTensor<2,2> get_stress_minus(const double,const double,SymmetricTensor<2,2>&);
template SymmetricTensor<2,3> get_stress_minus(const double,const double,SymmetricTensor<2,3>&);

