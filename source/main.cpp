#include "../include/parameter.h"
#include "../include/PhasefieldSMP.h"

int main(int argc, char *argv[]){

    try
      {
       using namespace dealii;

       if(argc==1)
       {
         std::cout<<"Please provide parameter file"<<std::endl;
         exit(0);
       }
       
       const std::string filename = argv[1];

       thesis::parameters::AllParameters param(filename);
       thesis::Phasefield<2> phasefield_2d(param);
       phasefield_2d.run(param);
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
