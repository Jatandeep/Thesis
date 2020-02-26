#include <iostream>
#include <string>
#include <fstream>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/tensor.h>
#include "../include/ElasticProblem.h"
#include "../include/parameter.h"
/*
#include <deal.II/base/symmetric_tensor.h>

template <int dim>
dealii::SymmetricTensor<2, 6,double> S_Tensor()
{

  dealii::SymmetricTensor<2, 6,double> tmp;
  for (unsigned int i = 0; i < 6; ++i)
    for (unsigned int j = 0; j < 6; ++j)
          tmp[i][j] = ((i == j) ? 1:0);
  return tmp;
}

template <int dim>
void Eigenvalue(){
   // dealii::SymmetricTensor<2,dim> unit_symmetric_tensor<dim>();
    // working

    const dealii::SymmetricTensor<2,6,double> s_tensor =
      S_Tensor<dim>();

    dealii::SymmetricTensor<2,dim,double> E;

    std::array <std::pair< double, dealii::Tensor< 1, 6, double > >,std::integral_constant< int, 6 >::value> eigen1;
    eigen1 = dealii::eigenvectors(s_tensor,dealii::SymmetricTensorEigenvectorMethod::ql_implicit_shifts);

    std::array <std::pair< double, dealii::Tensor< 1, dim, double > >,std::integral_constant< int, dim >::value> eigen2;
    eigen2 = dealii::eigenvectors(dealii::unit_symmetric_tensor<dim>(),dealii::SymmetricTensorEigenvectorMethod::ql_implicit_shifts);

    dealii::SymmetricTensor<4,dim,double> I;
    I=dealii::identity_tensor<dim>();

    double lambda =2;
    double mu =2;

    dealii::SymmetricTensor<4,dim,double> C_1;
    dealii::SymmetricTensor<4,dim,double> C_2;
    //C_2 = (dealii::unit_symmetric_tensor<dim>()*dealii::unit_symmetric_tensor<dim>());

    //std::array <std::pair< double, dealii::Tensor< 2, dim, double > >,std::integral_constant< int, dim >::value> eigens;
    //eigens = dealii::eigenvectors(dealii::identity_tensor<dim>(),dealii::SymmetricTensorEigenvectorMethod::ql_implicit_shifts);

    //C_2 = lambda*dealii::identity_tensor<dim>();            //values are zero?
    for (unsigned int i=0;i<dim;++i) {
        for (unsigned int j=0;j<dim;++j) {
            std::cout<<lambda*eigen2[i].second * eigen2[i].second * eigen2[j].second * eigen2[j].second<<std::endl;
            C_2[i][i][j][j] =  lambda*eigen2[i].second * eigen2[i].second * eigen2[j].second * eigen2[j].second;
            std::cout<<C_2[i][i][j][j]<<std::endl;
        }
    }

    for (unsigned int i=0;i<dim;++i) {
        for (unsigned int j=0;j<dim;++j) {
            //E[i][j]=eigen1[i].second;

        }

    }


    for (unsigned int i=0;i<dim;++i) {
        for (unsigned int j=0;j<dim;++j) {
            for (unsigned int m=0;m<dim;++m) {
                for (unsigned int n=0;n<dim;++n) {

                    //C_1[i][j][m][n] =E[i][j] * E[m][n];

                }

            }

        }

    }

}
*/

int main(int argc, char *argv[]){

    try
      {
       using namespace dealii;

       const std::string filename = argv[1];

       thesis::parameters::AllParameters param(filename);
       thesis::ElasticProblem<2> elastic_problem_2d(param);
       elastic_problem_2d.run(param);
         //Eigenvalue<3>();
      }
    catch (std::exception &exc)
      {
        std::cerr << std::endl << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        std::cerr << "Exception on processing: " << std::endl
                  << exc.what() << std::endl
                  << "Aborting!" << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        return 1;
      }
    catch (...)
      {
        std::cerr << std::endl << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        std::cerr << "Unknown exception!" << std::endl
                  << "Aborting!" << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        return 1;
      }
    return 0;
}
