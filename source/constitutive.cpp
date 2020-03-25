#include "../include/ElasticProblem.h"

using namespace dealii;

/*Gives BigC without spectral decomposition as in LKM lecture notes*/
template <int dim>
SymmetricTensor<4, dim> get_stress_strain_tensor(const double lambda
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

/*Gives stress at particular quadrature point*/
template<int dim>
SymmetricTensor<2,dim> get_stress(const double lambda
				 ,const double mu
				 ,SymmetricTensor<2,dim> &eps){

	SymmetricTensor<2,dim> stress;
	double tr_eps = trace(eps);

	stress = lambda*tr_eps*unit_symmetric_tensor<dim>() + 2*mu*eps;

return stress; 
}

/*Gives coefficient of first part of spectral BigC*/
double delsig_dellmbda(const double lambda
                      ,const double mu
		      ,unsigned int i
		      ,unsigned int j
		      ,unsigned int k
		      ,unsigned int l)
{
	double result = 0;

	result =  ((((i == k) && (j == l)) ? mu:0)
                   + (((i == l) && (j == k)) ? mu:0)
                   +((i==j) && (k==l)? lambda:0));
return result;
}

/*Gives coefficient of first part of positive BigC*/
template<int dim>
double delsig_plus_dellmbda(const double lambda
                           ,const double mu
                     	   ,unsigned int i
                     	   ,unsigned int j
			   ,unsigned int k
			   ,unsigned int l
                      	   ,SymmetricTensor<2,dim> & eps
			   ,double eps_eigenvalue)
{
        double result;
	double sgn_tr_eps, sgn_eps_eigenvalue;

	sgn_tr_eps = (trace(eps)>0) ? 1:0 ;
	sgn_eps_eigenvalue = (eps_eigenvalue>0) ? 1:0 ;

	result =  ((((i == k) && (j == l)) ? 0.5*mu*(1 + sgn_eps_eigenvalue):0)
		+ (((i == l) && (j == k)) ? 0.5*mu*(1 + sgn_eps_eigenvalue):0)
		+(((i==j) && (k==l)) ? (0.5*lambda*(1 + sgn_tr_eps)):0));
return result;
}

/*Gives coefficient of first part of negative BigC*/
template<int dim>
double delsig_minus_dellmbda(const double lambda
	                    ,const double mu
                     	    ,unsigned int i
                      	    ,unsigned int j
			    ,unsigned int k
			    ,unsigned int l
                            ,SymmetricTensor<2,dim> & eps
                            ,double eps_eigenvalue)
{
        double result;
        double sgn_tr_eps, sgn_eps_eigenvalue;

        sgn_tr_eps = (trace(eps)>0) ? 1:0 ;
        sgn_eps_eigenvalue = (eps_eigenvalue>0) ? 1:0 ;

	result =  ((((i == k) && (j == l)) ? 0.5*mu*(1 - sgn_eps_eigenvalue):0)
		+(((i == l) && (j == k)) ? 0.5*mu*(1 - sgn_eps_eigenvalue):0)
		+(((i==j) && (k==l))? (0.5*lambda*(1 - sgn_tr_eps)):0));
	
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
double get_stress_plus_eigenvalue(const double lambda
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
double get_stress_minus_eigenvalue(const double lambda
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

	SymmetricTensor<4,dim> C_1_plus;
       	SymmetricTensor<4,dim> C_1_minus;
	SymmetricTensor<4,dim> C_2_plus;
        SymmetricTensor<4,dim> C_2_minus;

	SymmetricTensor<4,dim> C_1;
	SymmetricTensor<4,dim> C_2;

	std::array <std::pair< double, Tensor< 1, dim, double > >,std::integral_constant< int, dim >::value> eigen;
  	eigen = eigenvectors(eps,SymmetricTensorEigenvectorMethod::ql_implicit_shifts);


        for (unsigned int i=0;i<dim;++i) 
                for (unsigned int j=i;j<dim;++j)
			for(unsigned int k=0;k<dim;++k )
				for(unsigned int l=k;l<dim;++l) {
                        C_1_plus[i][j][k][l] =  delsig_plus_dellmbda(lambda,mu,i,j,k,l,eps,eigen[i].first) * eigen[i].second * eigen[i].second * eigen[k].second * eigen[k].second;

                }
        

        for (unsigned int i=0;i<dim;++i) 
                for (unsigned int j=i;j<dim;++j)
			for(unsigned int k=0;k<dim;++k)
				for(unsigned int l=k;l<dim;++l) {
                        C_1_minus[i][j][k][l] =  delsig_minus_dellmbda(lambda,mu,i,j,k,l,eps,eigen[i].first) * eigen[i].second * eigen[i].second * eigen[k].second * eigen[k].second;

                }
       

        for (unsigned int i=0;i<dim;++i) {
                for (unsigned int j=0;j<dim;++j) {
                if( !(eigen[i].first==eigen[j].first) ){
                        C_2_plus[i][j][i][j] =0.5* ( ( get_stress_plus_eigenvalue(lambda,mu,eps,eigen[j].first) - get_stress_plus_eigenvalue(lambda,mu,eps,eigen[i].first))/(eigen[j].first - eigen[i].first) ) *( eigen[i].second * eigen[j].second * eigen[i].second * eigen[j].second + eigen[i].second * eigen[j].second * eigen[j].second * eigen[i].second );

                }
                }
        }

        for (unsigned int i=0;i<dim;++i) {
                for (unsigned int j=0;j<dim;++j) {
                if( !(eigen[i].first==eigen[j].first) ){
                        C_2_minus[i][j][i][j] =0.5* ( ( get_stress_minus_eigenvalue(lambda,mu,eps,eigen[j].first) - get_stress_minus_eigenvalue(lambda,mu,eps,eigen[i].first))/(eigen[j].first - eigen[i].first) ) *( eigen[i].second * eigen[j].second * eigen[i].second * eigen[j].second +eigen[i].second * eigen[j].second * eigen[j].second * eigen[i].second );

                }
                }
        }




  	for (unsigned int i = 0;i < dim; ++i) 
		for (unsigned int j = i;j < dim; ++j)
		     for (unsigned int k = 0; k < dim; ++k)
       			for (unsigned int l = k; l < dim; ++l)			 {
          		C_1[i][j][k][l] =  delsig_dellmbda(lambda,mu,i,j,k,l) * eigen[i].second * eigen[i].second * eigen[k].second * eigen[k].second;

      		}
 

        for (unsigned int i=0;i<dim;++i) {
                for (unsigned int j=0;j<dim;++j) {
		if( !(eigen[i].first==eigen[j].first) ){
                        C_2[i][j][i][j] =0.5* ( ( get_stress_eigenvalue(lambda,mu,eps,eigen[j].first) - get_stress_eigenvalue(lambda,mu,eps,eigen[i].first))/(eigen[j].first - eigen[i].first) ) *( eigen[i].second * eigen[j].second * eigen[i].second * eigen[j].second +eigen[i].second * eigen[j].second * eigen[j].second * eigen[i].second );

                }
		}
        }

/*COMMENT- C_2  or (C_2_plus+C_2_minus) part doesn't contribute to the result*/

//return(C_1+C_2);
return (C_1_plus + C_1_minus + C_2_plus + C_2_minus);
}


template SymmetricTensor<4,2> get_stress_strain_tensor(double,double);
template SymmetricTensor<4,3> get_stress_strain_tensor(double,double);
template SymmetricTensor<4,2> get_BigC(double,double,SymmetricTensor<2,2>&);
template SymmetricTensor<4,3> get_BigC(double,double,SymmetricTensor<2,3>&);
template SymmetricTensor<2,2> get_stress(const double,const double,SymmetricTensor<2,2>&);
template SymmetricTensor<2,3> get_stress(const double,const double,SymmetricTensor<2,3>&);
