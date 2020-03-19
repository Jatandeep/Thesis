#include "../include/ElasticProblem.h"

using namespace dealii;


template <int dim>
SymmetricTensor<4, dim> get_stress_strain_tensor(const double lambda,
                                                 const double mu)
{

  SymmetricTensor<4, dim> C;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        for (unsigned int l = k; l < dim; ++l)
          C[i][j][k][l] = ((((i == k) && (j == l)) ? mu:0)
                            + (((i == l) && (j == k)) ? mu:0)
                             +((i==j) && (k==l)? lambda:0));


  return (C);
}


template<int dim>
double delsig_dellmbda(const double lambda
                      ,const double mu
		      ,unsigned int i
		      , unsigned int j
		      ,SymmetricTensor<2,dim> &epsilon)
{
	double value;

	if(i==j)
	value = lambda + 2*mu;
	else
	value = lambda;

return value;
}

template<int dim>
double get_stress_eigenvalue(const double lambda
                      ,const double mu
                      ,SymmetricTensor<2,dim> &epsilon
			,double strain_eigenvalue)
{
	double value;
	value = lambda*trace(epsilon) + 2*mu*strain_eigenvalue;

return value;
}


template <int dim>
SymmetricTensor<4,dim> get_BigC(const double lambda
                                        ,const double mu
                                        ,SymmetricTensor<2,dim> &epsilon)
{

	SymmetricTensor<4,dim> C_1;
	SymmetricTensor<4,dim> C_2;

	std::array <std::pair< double, dealii::Tensor< 1, dim, double > >,std::integral_constant< int, dim >::value> eigen;
  	eigen = dealii::eigenvectors(epsilon/*dealii::unit_symmetric_tensor<dim>()*/,
				     dealii::SymmetricTensorEigenvectorMethod::ql_implicit_shifts);

  	for (unsigned int i=0;i<dim;++i) {
      		for (unsigned int j=0;j<dim;++j) {
          		C_1[i][i][j][j] =  delsig_dellmbda(lambda,mu,i,j,epsilon) * eigen[i].second * eigen[i].second * eigen[j].second * eigen[j].second;

      		}
  	}

        for (unsigned int i=0;i<dim;++i) {
                for (unsigned int j=0;j<dim;++j) {
		if(!(i==j) && !(eigen[i].first==eigen[j].first) ){
                        C_2[i][i][j][j] =0.5* ( ( get_stress_eigenvalue(lambda,mu,epsilon,eigen[j].first) - get_stress_eigenvalue(lambda,mu,epsilon,eigen[i].first))/(eigen[j].first - eigen[i].first) ) *( eigen[i].second * eigen[j].second * eigen[i].second * eigen[j].second +  eigen[i].second * eigen[j].second * eigen[j].second * eigen[i].second );

                }
		}
        }

return (C_1+C_2);

}






template SymmetricTensor<4,2> get_stress_strain_tensor(double,double);
template SymmetricTensor<4,3> get_stress_strain_tensor(double,double);
template SymmetricTensor<4,2> get_BigC(double,double,SymmetricTensor<2,2>&);
template SymmetricTensor<4,3> get_BigC(double,double,SymmetricTensor<2,3>&);
