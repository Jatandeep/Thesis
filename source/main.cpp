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

       /*Extracting filename for outputs*/
       std::size_t found = filename.find_last_of("/\\");
       std::string rawfilename = filename.substr(found+1).substr(0, filename.substr(found+1).find_last_of("."));
       
       thesis::parameters::AllParameters param(filename);
       thesis::Phasefield<2> phasefield_2d(param);
       phasefield_2d.run(param,rawfilename);
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
