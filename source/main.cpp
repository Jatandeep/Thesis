#include <iostream>
#include <string>
#include <fstream>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/tensor.h>
#include "../include/ElasticProblem.h"
#include "../include/parameter.h"





int main(int argc, char *argv[]){

    try
      {
        using namespace dealii;

        const std::string filename = argv[1];

      thesis::AllParameters parameter(filename);
      thesis::ElasticProblem<2> elastic_problem_2d(parameter);
//        thesis::ElasticProblem<2> elastic_problem_2d(filename);
        elastic_problem_2d.run(parameter);

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
