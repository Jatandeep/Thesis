#include <iostream>
#include <string>
#include <fstream>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/tensor.h>
#include "../include/parameter.h"
#include "../include/Phasefield.h"
//#include "../include/PhasefieldSMP.h"
//#include "../include/ElasticProblem.h"

int main(int argc, char *argv[]){

    try
      {
       using namespace dealii;
	
       const std::string filename = argv[1];

       thesis::parameters::AllParameters param(filename);
       thesis::Phasefield<2> phasefield_2d(param);
       phasefield_2d.run(param);
//	thesis::ElasticProblem<2> elastic_2d(param);
//	elastic_2d.run(param);
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
