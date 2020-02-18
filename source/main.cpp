#include <iostream>
#include <string>
#include <fstream>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/tensor.h>
#include "../include/lin_elas.h"
#include "../include/parameter.h"
#include <deal.II/base/timer.h>


template <int dim>
void ElasticProblem<dim>::run(){

    using namespace dealii;

    Timer timer;
    timer.start();

    for (unsigned int cy =0; cy<parameter.cycles; ++cy)
      {
        std::cout << "Cycle " << cy << ':' << std::endl;
        if (cy == 0)
          {
            //GridGenerator::hyper_cube (triangulation, -1, 1);
            import_mesh();
            triangulation.refine_global (parameter.gl_ref);
          }
        else
          refine_grid ();
        std::cout << "   Number of active cells:       "
                  << triangulation.n_active_cells()
                  << std::endl;
        setup_system ();
        std::cout << "   Number of degrees of freedom: "
                  << dof_handler.n_dofs()
                  << std::endl;
        assemble_system ();
        solve ();
        output_results (cy);
        std::cout<<"time:"<<timer()<<" solution.norm():"<<solution.l2_norm()<<std::endl;
      }
}

int main(int argc, char *argv[]){

    try
      {
        using namespace dealii;

        const std::string filename = argv[1];
        ElasticProblem<2> elastic_problem_2d(filename);
        elastic_problem_2d.run();


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
