#include "../include/utilities.h"

using namespace dealii;

/*To get tensile/positive component of stress to be used in assembly of u*/
/*Gives stress_plus at particular quadrature point*//*Thesis_report:Equation 2.38*/
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

return stress_plus; 
}

/*To get compressive/negative component of stress to be used in assembly of u*/
/*Gives stress_minus at particular quadrature point*//*Thesis_report:Equation 2.39*/
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

return stress_minus; 
}

//Gives an array of +ve stress eigenvalues/*To be used in get_BigC_plus()*/
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

//Gives an array of -ve stress eigenvalues/*To be used in get_BigC_minus()*/
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

//Gives part of +ve coefficient after applying L'hospital rule (in second part of +ve BigC)/*To be used in get_BigC_plus()*/
template<int dim>
double delsig_dellmbda_b_plus(const double lambda
                             ,const double mu
                     	     ,SymmetricTensor<2,dim> & eps
			     ,double eps_eigenvalue)
{
    double result;
	double sgn_tr_eps, sgn_eps_eigenvalue;

	sgn_tr_eps = get_sign(trace(eps));

	sgn_eps_eigenvalue = get_sign(eps_eigenvalue);

	result = 0.5*lambda*(1+sgn_tr_eps) + mu*(1+sgn_eps_eigenvalue);

return result;
}

//Gives part of -ve coefficient after applying L'hospital rule (in second part of -ve BigC)/*To be used in get_BigC_minus()*/
template<int dim>
double delsig_dellmbda_b_minus(const double lambda
                             ,const double mu
                     	     ,SymmetricTensor<2,dim> & eps
			     ,double eps_eigenvalue)
{
    double result;
	double sgn_tr_eps, sgn_eps_eigenvalue;

	sgn_tr_eps = get_sign(trace(eps));

	sgn_eps_eigenvalue = get_sign(eps_eigenvalue);
	
	result = 0.5*lambda*(1-sgn_tr_eps) + mu*(1-sgn_eps_eigenvalue);

return result;
}

/*Gives coefficient of first part of positive BigC*//*To be used in get_BigC_plus()*/
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

	sgn_tr_eps = get_sign(trace(eps));	

	sgn_eps_eigenvalue = get_sign(eps_eigenvalue);

	result = 0.5*lambda*(1+sgn_tr_eps) + ((i==j) ? mu*(1+sgn_eps_eigenvalue):0);

return result;
}

/*Gives coefficient of first part of negative BigC*//*To be used in get_BigC_minus()*/
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

	sgn_tr_eps = get_sign(trace(eps));
	sgn_eps_eigenvalue = get_sign(eps_eigenvalue);

      
	result = 0.5*lambda*(1-sgn_tr_eps) + ((i==j) ? mu*(1-sgn_eps_eigenvalue):0);

return result;
}

/*To get positive component of BigC to be used in assembly of u*//*Thesis_report:Equation 3.43*/
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

return C_total_plus;
}

/*To get negative component of BigC to be used in assembly of u*//*Thesis_report:Equation 3.48*/
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

return C_total_minus;
}


template SymmetricTensor<4,2> get_BigC_plus(double,double,SymmetricTensor<2,2>&);
template SymmetricTensor<4,3> get_BigC_plus(double,double,SymmetricTensor<2,3>&);
template SymmetricTensor<4,2> get_BigC_minus(double,double,SymmetricTensor<2,2>&);
template SymmetricTensor<4,3> get_BigC_minus(double,double,SymmetricTensor<2,3>&);
template SymmetricTensor<2,2> get_stress_plus(const double,const double,SymmetricTensor<2,2>&);
template SymmetricTensor<2,3> get_stress_plus(const double,const double,SymmetricTensor<2,3>&);
template SymmetricTensor<2,2> get_stress_minus(const double,const double,SymmetricTensor<2,2>&);
template SymmetricTensor<2,3> get_stress_minus(const double,const double,SymmetricTensor<2,3>&);

