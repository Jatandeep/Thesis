#include"../include/Phasefield.h"

using namespace dealii;
using namespace thesis;
using namespace parameters;

template <int dim>
void Phasefield<dim>::solve_linear_steps(const parameters::AllParameters &param
                             ,dealii::BlockVector<double> &solution)
{

  implement_constraints(param,current_time_m);
	
  assemble_d(param,solution_m);
  solve_system(param,solution_m);

  assemble_u(param,solution_m);
  solve_system(param,solution_m);


}

