#include "../include/ElasticProblem.h"
#include "../include/others.h"
#include "../include/constants.h"
#include "../include/constitutive.h"
#include "../include/utilities.h"

using namespace dealii;
using namespace thesis;
using namespace parameters;

template <int dim>
ElasticProblem<dim>::ElasticProblem(const AllParameters &param)
  : fe_m(FE_Q<dim>(param.fesys.fe_degree),dim)
  , quadrature_formula_m(param.fesys.fe_degree+1)
  , disp_extractor(0)
{
    import_mesh(param);
    dof_handler_m.initialize(triangulation_m, fe_m);
    setup_system();
}



template <int dim>
ElasticProblem<dim>::~ElasticProblem ()
{
  dof_handler_m.clear ();
}


template <int dim>
void ElasticProblem<dim>::setup_system ()
{

  dof_handler_m.distribute_dofs (fe_m);
  //DoFRenumbering::Cuthill_McKee(dof_handler_m);
  constraints_m.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler_m,
                                           constraints_m);
  constraints_m.close ();
  tangent_matrix_m.clear();

  DynamicSparsityPattern dsp(dof_handler_m.n_dofs(), dof_handler_m.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler_m
                                  ,dsp
                                  ,constraints_m
                                  ,true);/*keep_constrained_dofs = */
  sparsity_pattern_m.copy_from (dsp);
  tangent_matrix_m.reinit (sparsity_pattern_m);
  system_rhs_m.reinit (dof_handler_m.n_dofs());
  solution_m.reinit (dof_handler_m.n_dofs());

  std::cout << "   Number of active cells:       "
                            << triangulation_m.n_active_cells()
                            << std::endl;
  std::cout << "   Number of degrees of freedom: "
                            << dof_handler_m.n_dofs()
                            << std::endl;
}

template <int dim>
void ElasticProblem<dim>::assemble_system (const AllParameters &param,Vector<double> & newton_update)

{

  FEValues<dim> fe_values (fe_m, quadrature_formula_m,
                           update_values   | update_gradients |
                           update_quadrature_points | update_JxW_values);
  const unsigned int   dofs_per_cell = fe_m.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula_m.size();
  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  SymmetricTensor<4,dim> BigC;

  SymmetricTensor<4,dim> Big_C;

  Big_C = get_const_BigC<dim>(param.materialmodel.lambda,param.materialmodel.mu);

  for (const auto &cell : dof_handler_m.active_cell_iterators())
    {
      cell_matrix = 0;
      cell_rhs=0;

      fe_values.reinit(cell);

      cell->get_dof_indices (local_dof_indices);

      std::vector<SymmetricTensor<2,dim>> epsilon_vals(n_q_points);
      fe_values[disp_extractor].get_function_symmetric_gradients(newton_update,epsilon_vals);

        for (unsigned int q = 0; q < n_q_points; ++q){

	    BigC = get_BigC(param.materialmodel.lambda,param.materialmodel.mu,epsilon_vals[q]);
 
	    SymmetricTensor<2,dim> sigma = get_stress(param.materialmodel.lambda
						     ,param.materialmodel.mu
						     ,epsilon_vals[q]);

            for (unsigned int i = 0; i < dofs_per_cell; ++i){

                const SymmetricTensor<2, dim> sym_grad_shape_i = fe_values[disp_extractor].symmetric_gradient(i, q);

                cell_rhs(i) -= (sigma * sym_grad_shape_i)*fe_values.JxW(q);

                for (unsigned int j = i; j < dofs_per_cell; ++j){
                    const SymmetricTensor<2, dim> sym_grad_shape_j = fe_values[disp_extractor].symmetric_gradient(j, q);

                    cell_matrix(i, j) += (sym_grad_shape_i *
                                          BigC *
                                          sym_grad_shape_j
                                          ) *
                                         fe_values.JxW(q);

                }
            }
          }

        for(unsigned int i=0;i<dofs_per_cell;++i)
            for(unsigned int j=0;j<i;++j)
                cell_matrix(i,j)=cell_matrix(j,i);

       constraints_m.distribute_local_to_global(cell_matrix,cell_rhs,
                                  local_dof_indices,tangent_matrix_m,system_rhs_m,false);  //copy local to global
  }


}

template <int dim>
void ElasticProblem<dim>::assemble_body_forces(const AllParameters &param)

{

  FEValues<dim> fe_values (fe_m, quadrature_formula_m,
                           update_values   | update_gradients |
                           update_quadrature_points | update_JxW_values);
  const unsigned int   dofs_per_cell = fe_m.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula_m.size();
  Vector<double>       cell_rhs (dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  std::vector<Vector<double> >      rhs_values (n_q_points,Vector<double>(dim+1));
  const Others<dim> right_hand_side;

  for (const auto &cell : dof_handler_m.active_cell_iterators())
    {

      cell_rhs = 0;

      fe_values.reinit(cell);

      double step_fraction = double(current_time)/double(param.time.end_time);
      right_hand_side.vector_value_list(fe_values.get_quadrature_points(), rhs_values);

         for (unsigned int q = 0; q < n_q_points; ++q)
            for (unsigned int i = 0; i < dofs_per_cell; ++i){

                const unsigned int
                component_i = fe_m.system_to_component_index(i).first;
                cell_rhs(i) += fe_values.shape_value(i,q) *
                               rhs_values[q][component_i] *step_fraction*
                               fe_values.JxW(q);
            }

        cell->get_dof_indices (local_dof_indices);

    constraints_m.distribute_local_to_global(cell_rhs,
                                  local_dof_indices,
                                  system_rhs_m);//copy local to global
  }

}

template <int dim>
void ElasticProblem<dim>::solve_nonlinear_newton(const AllParameters &param,
                                                 Vector<double> &solution_delta){

    Vector<double> newton_update(dof_handler_m.n_dofs());

    error_residual.reset();
    error_residual_0.reset();
    error_residual_norm.reset();

    print_header();
    unsigned int new_iter = 0;
    for (; new_iter < param.newtonraphson.max_new_ite; ++new_iter) {

        std::cout << " " << std::setw(2) << new_iter << " " << std::flush;

        tangent_matrix_m = 0.0;
        system_rhs_m = 0.0;

        make_constraints(new_iter);

        assemble_system(param,newton_update);
        assemble_body_forces(param);

        get_error_residual(error_residual);

        if(new_iter==0)
            error_residual_0 = error_residual;

        error_residual_norm = error_residual;
        error_residual_norm.normalize(error_residual_0);

        if(new_iter > 0 && error_residual_norm.u < param.newtonraphson.res_tol){
            std::cout<<"Converged"<<std::endl;
            print_footer();
            break;
        }

        const std::pair<unsigned int,double>
                lin_solver_output = solve_linear_sys(param,newton_update);

        solution_delta += newton_update;

        std::cout << " | " << std::fixed << std::setprecision(3) << std::setw(7)
                            << std::scientific << lin_solver_output.first << "  "
                            << lin_solver_output.second << "  " << error_residual_norm.u
                            << "  " << std::endl;

    }
    AssertThrow (new_iter < param.newtonraphson.max_new_ite,
                   ExcMessage("No convergence in nonlinear solver!"));

}

template <int dim>
std::pair<unsigned int,double> ElasticProblem<dim>::solve_linear_sys(const AllParameters &param,Vector<double> &newton_update)
{

  unsigned int lin_ite = 0;
  double lin_res = 0;
  newton_update = 0;    // (solution)

  SolverControl           solver_control (tangent_matrix_m.m(),
                                          param.linearsolver.cg_tol*system_rhs_m.l2_norm());
  GrowingVectorMemory<Vector<double> > GVM;
  SolverCG<Vector<double>>    cg (solver_control,GVM);
  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(tangent_matrix_m, param.linearsolver.relax_prm);
  cg.solve (tangent_matrix_m, newton_update, system_rhs_m,
            preconditioner);

  lin_ite = solver_control.last_step();
  lin_res = solver_control.last_value();

  constraints_m.distribute (newton_update);
  return std::make_pair(lin_ite,lin_res);
}

template <int dim>
void ElasticProblem<dim>::refine_grid (const AllParameters &param)
{
  Vector<float> estimated_error_per_cell (triangulation_m.n_active_cells());
  KellyErrorEstimator<dim>::estimate (dof_handler_m,
                                      QGauss<dim-1>(param.fesys.quad_order),
                                      typename FunctionMap<dim>::type(),
                                      solution_m,
                                      estimated_error_per_cell);
  GridRefinement::refine_and_coarsen_fixed_number (triangulation_m,
                                                   estimated_error_per_cell,
                                                   param.geometrymodel.act_ref, param.geometrymodel.act_cors);
  triangulation_m.execute_coarsening_and_refinement ();

}


template <int dim>
void ElasticProblem<dim>::output_results (const double cycle) const
{

  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler_m);
  std::vector<std::string> solution_names;
  switch (dim)
    {
    case 1:
      solution_names.emplace_back("displacement");
      break;
    case 2:
      solution_names.emplace_back("x_displacement");
      solution_names.emplace_back("y_displacement");
      break;
    case 3:
      solution_names.emplace_back("x_displacement");
      solution_names.emplace_back("y_displacement");
      solution_names.emplace_back("z_displacement");
      break;
    default:
      Assert (false, ExcNotImplemented());
    }
  data_out.add_data_vector (solution_m, solution_names);
  data_out.build_patches ();
  std::ofstream output ("solution-" + std::to_string(cycle) + ".vtk");
  data_out.write_vtk (output);
}

template <int dim>
void ElasticProblem<dim>::import_mesh(const AllParameters &param){

    std::string grid_name;
    grid_name += param.geometrymodel.meshfile;

    GridIn<dim> grid_in;
    grid_in.attach_triangulation(triangulation_m);
    std::ifstream input_file(grid_name.c_str());
    Assert(dim==2, ExcInternalError());
    grid_in.read_abaqus(input_file,false);

    triangulation_m.refine_global (param.geometrymodel.gl_ref);
}



template <int dim>
void ElasticProblem<dim>::print_header(){
    const unsigned int l_width = 80;
        for (unsigned int i = 0; i < l_width; ++i)
        {
            std::cout << "_";
        }
        std::cout << std::endl;

        std::cout << "           SOLVER STEP            "
                    << " |  LIN_IT   LIN_RES    RES_NORM    "
                    << std::endl;

        for (unsigned int i = 0; i < l_width; ++i)
        {
            std::cout << "_";
        }
        std::cout << std::endl;

}

template <int dim>
void ElasticProblem<dim>::print_footer(){
    const unsigned int l_width = 80;
        for (unsigned int i = 0; i < l_width; ++i)
        {
            std::cout << "_";
        }
        std::cout << std::endl;

        std::cout << "Errors:" << std::endl
                    << "Rhs: \t\t" << error_residual_norm.u << std::endl
                    << std::endl;

}

template <int dim>
void ElasticProblem<dim>::get_error_residual(Error& error_residual){


    Vector<double> err_res(dof_handler_m.n_dofs());
    for (unsigned int i=0;i<dof_handler_m.n_dofs();++i) {
        if(!constraints_m.is_constrained(i)){
            err_res(i)=system_rhs_m(i);
           }

        }
    error_residual.u = err_res.l2_norm();

}

template <int dim>
void ElasticProblem<dim>::make_constraints(unsigned int &itr){
    std::cout<<" CST "<<std::flush;

    if(itr>=1){
        return;
    }

   constraints_m.clear();

   VectorTools::interpolate_boundary_values (dof_handler_m,
                                              0,
                                              Functions::ZeroFunction<dim>(dim),
                                              constraints_m
                                              ,fe_m.component_mask(disp_extractor));

    constraints_m.close ();
}


template <int dim>
void ElasticProblem<dim>::run(const AllParameters &param){


    using namespace constants;

    Vector<double>       solution_delta(dof_handler_m.n_dofs());

    long double present_time_tol;
    present_time_tol = time_tol_fac * param.time.delta_t;
    current_time = param.time.start_time;

    while (current_time < param.time.end_time + present_time_tol)
    {
        if(current_time > param.time.delta_t + present_time_tol)
/*	To be implemented
	{
            	refine_grid(param);
		setup_system();
		solution_delta.reinit(dof_handler_m.n_dofs());
	}
*/
        solution_delta = 0.0;
        solve_nonlinear_newton(param,solution_delta);
        solution_m = solution_delta;
        output_results(current_time);
	
        std::cout<<" solution.norm():"<<solution_m.l2_norm()<<std::endl;
//	std::cout<<"solution at node(1000):"<<solution_m(1000)<<std::endl;
        current_time += param.time.delta_t;
    }

}

template class thesis::ElasticProblem<2>;
template class thesis::ElasticProblem<3>;

