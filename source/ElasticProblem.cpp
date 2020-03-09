#include "../include/ElasticProblem.h"
#include "../include/others.h"
#include "../include/HyperCubeWithRefinedHole.h"

using namespace dealii;
using namespace thesis;
using namespace parameters;

template <int dim>
SymmetricTensor<4, dim> get_stress_strain_tensor(const double lambda,
                                                 const double mu)
{

  SymmetricTensor<4, dim> C_1;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        for (unsigned int l = k; l < dim; ++l)
          C_1[i][j][k][l] = ((((i == k) && (j == l)) ? mu:0)
                            + (((i == l) && (j == k)) ? mu:0)
                             +((i==j) && (k==l)? lambda:0));

/*
  SymmetricTensor<4,dim,double> C_2;
  std::array <std::pair< double, dealii::Tensor< 1, dim, double > >,std::integral_constant< int, dim >::value> eigen2;
  eigen2 = dealii::eigenvectors(dealii::unit_symmetric_tensor<dim>(),dealii::SymmetricTensorEigenvectorMethod::ql_implicit_shifts);

  for (unsigned int i=0;i<dim;++i) {
      for (unsigned int j=0;j<dim;++j) {
          C_2[i][i][j][j] =  lambda * eigen2[i].second * eigen2[i].second * eigen2[j].second * eigen2[j].second;

      }
  }
*/
  return (C_1/*+C_2*/);
}

template <int dim>
const SymmetricTensor<4, dim> ElasticProblem<dim>::stress_strain_tensor =
  get_stress_strain_tensor<dim>(1, 1);

template <int dim>
ElasticProblem<dim>::ElasticProblem(AllParameters param):
    fe(FE_Q<dim>(param.fesys.fe_degree),dim), dof_handler(triangulation),
    quadrature_formula(param.fesys.fe_degree+1)
{}



template <int dim>
ElasticProblem<dim>::~ElasticProblem ()
{
  dof_handler.clear ();
}


template <int dim>
void ElasticProblem<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);
  //DoFRenumbering::Cuthill_McKee(dof_handler);
  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
                                           hanging_node_constraints);
  hanging_node_constraints.close ();
  tangent_matrix.clear();
  DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler,
                                  dsp,
                                  hanging_node_constraints,
                                  /*keep_constrained_dofs = */ true);
  sparsity_pattern.copy_from (dsp);
  tangent_matrix.reinit (sparsity_pattern);
  system_rhs.reinit (dof_handler.n_dofs());
  solution_delta.reinit(dof_handler.n_dofs());
  solution.reinit (dof_handler.n_dofs());

}

template <int dim>
void ElasticProblem<dim>::assemble_system_matrix (Error & error_residual,
                                                  dealii::Vector<double> & newton_update,
                                                  unsigned int &new_ite,bool error)
{

  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values   | update_gradients |
                           update_quadrature_points | update_JxW_values);
  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();
  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  const FEValuesExtractors::Vector disp(0);


  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell_matrix = 0;

      fe_values.reinit(cell);

        for (unsigned int q = 0; q < n_q_points; ++q)
            for (unsigned int i = 0; i < dofs_per_cell; ++i){
                const SymmetricTensor<2, dim> phi_i_u = fe_values[disp].symmetric_gradient(i, q);

                for (unsigned int j = 0; j < dofs_per_cell; ++j){
                    const SymmetricTensor<2, dim> phi_j_u = fe_values[disp].symmetric_gradient(j, q);

                    cell_matrix(i, j) += (phi_i_u *
                                          stress_strain_tensor *
                                          phi_j_u
                                          ) *
                                         fe_values.JxW(q);

                }
          }

    cell->get_dof_indices (local_dof_indices);

    hanging_node_constraints.distribute_local_to_global(cell_matrix,
                                  local_dof_indices,
                                  tangent_matrix);//copy local to global
  }

  if(error){

      FullMatrix<double> K(dof_handler.n_dofs(),dof_handler.n_dofs());
      Vector<double> u(dof_handler.n_dofs());

      for (unsigned int i=0;i<dof_handler.n_dofs();++i) {
          if(!hanging_node_constraints.is_constrained(i)){
              u(i)=newton_update(i);
             }
      }

      K.copy_from(tangent_matrix);
      Vector<double> tmp(dof_handler.n_dofs());
        K.vmult(tmp,u);

        if(new_ite==0)
            error_residual.u = tmp.l2_norm();
        else {
            error_residual.u -= tmp.l2_norm();
        }
        std::cout<<"error_residual_matrix:"<<error_residual.u<<std::endl;


  }
}

template <int dim>
void ElasticProblem<dim>::assemble_neumann_forces(Error & error_residual,AllParameters param
                                                  ,unsigned int &new_ite,bool error)
{

  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values   | update_gradients |
                           update_quadrature_points | update_JxW_values);
  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();
  Vector<double>       cell_rhs (dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  std::vector<Vector<double> >      rhs_values (n_q_points,Vector<double>(dim+1));
  const Others<dim> right_hand_side;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {

      cell_rhs = 0;

      fe_values.reinit(cell);

    double step_fraction = double(current_time_step)/double(param.fesys.n_time_steps);
      right_hand_side.vector_value_list(fe_values.get_quadrature_points(), rhs_values);


         for (unsigned int q = 0; q < n_q_points; ++q)
            for (unsigned int i = 0; i < dofs_per_cell; ++i){

                const unsigned int
                component_i = fe.system_to_component_index(i).first;
                cell_rhs(i) += fe_values.shape_value(i,q) *
                               rhs_values[q][component_i] *step_fraction*
                               fe_values.JxW(q);
            }

        cell->get_dof_indices (local_dof_indices);

    hanging_node_constraints.distribute_local_to_global(cell_rhs,
                                  local_dof_indices,
                                  system_rhs);//copy local to global
  }
  if(error){

      if(new_ite==0){

      Vector<double> err_res(dof_handler.n_dofs());
      for (unsigned int i=0;i<dof_handler.n_dofs();++i) {
          if(!hanging_node_constraints.is_constrained(i)){
              err_res(i)=system_rhs(i);
             }
          }

      error_residual.u += err_res.l2_norm();
      std::cout<<"error_residual_rhs:"<<error_residual.u<<std::endl;
      std::cout<<"system_rhs:"<<system_rhs.l2_norm()<<std::endl;
    }
  }
}

template <int dim>
void ElasticProblem<dim>::solve_nonlinear_newton(AllParameters param,
                                                 Vector<double> &solution_delta){

    Vector<double> newton_update(dof_handler.n_dofs());

    error_residual.reset();
    error_residual_0.reset();
    error_residual_norm.reset();

    print_header();
    unsigned int new_ite = 0;
    for (; new_ite < param.fesys.max_new_ite; ++new_ite) {

        std::cout << " " << std::setw(2) << new_ite << " " << std::flush;

        tangent_matrix = 0.0;
        system_rhs = 0.0;
        make_constraints(new_ite);

        bool error = false;
        //bool error = true;
        assemble_system_matrix(error_residual,newton_update,new_ite,error);
        assemble_neumann_forces(error_residual,param,new_ite,error);

        get_error_residual(error_residual,solution_delta);

        if(new_ite==0)
            error_residual_0 = error_residual;

        error_residual_norm = error_residual;
        error_residual_norm.normalize(error_residual_0);

        if(new_ite > 0 && error_residual_norm.u <= param.fesys.res_tol){
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
    AssertThrow (new_ite < param.fesys.max_new_ite,
                   ExcMessage("No convergence in nonlinear solver!"));

}

template <int dim>
std::pair<unsigned int,double> ElasticProblem<dim>::solve_linear_sys(AllParameters param,Vector<double> &newton_update)
{

  unsigned int lin_ite = 0;
  double lin_res = 0;
  newton_update = 0;    // (solution)

  SolverControl           solver_control (/*param.fesys.steps ,param.fesys.tol*/
                                          tangent_matrix.m(),param.fesys.tol*system_rhs.l2_norm());
  GrowingVectorMemory<Vector<double> > GVM;
  SolverCG<Vector<double>>    cg (solver_control,GVM);
  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(tangent_matrix, param.fesys.relax_prm);
  cg.solve (tangent_matrix, newton_update, system_rhs,
            preconditioner);

  lin_ite = solver_control.last_step();
  lin_res = solver_control.last_value();

  hanging_node_constraints.distribute (newton_update);
  return std::make_pair(lin_ite,lin_res);
}

template <int dim>
void ElasticProblem<dim>::refine_grid (AllParameters param)
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
  KellyErrorEstimator<dim>::estimate (dof_handler,
                                      QGauss<dim-1>(param.fesys.quad_order),
                                      typename FunctionMap<dim>::type(),
                                      solution,
                                      estimated_error_per_cell);
  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                   estimated_error_per_cell,
                                                   param.fesys.act_ref, param.fesys.act_cors);
  triangulation.execute_coarsening_and_refinement ();

}



template <int dim>
void ElasticProblem<dim>::output_results (const unsigned int cycle) const
{

  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
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
  data_out.add_data_vector (solution, solution_names);
  data_out.build_patches ();
  std::ofstream output ("solution-" + std::to_string(cycle) + ".vtk");
  data_out.write_vtk (output);
}

template <int dim>
void ElasticProblem<dim>::import_mesh(AllParameters param){

    std::string grid_name;
    grid_name += param.geometrymodel.meshfile;

    GridIn<dim> grid_in;
    grid_in.attach_triangulation(triangulation);
    std::ifstream input_file(grid_name.c_str());
    Assert(dim==2, ExcInternalError());
    grid_in.read_abaqus(input_file,false);

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
void ElasticProblem<dim>::get_error_residual(Error& error_residual,Vector<double> &newton_update){

/*
    Vector<double> err_res(dof_handler.n_dofs());
    for (unsigned int i=0;i<dof_handler.n_dofs();++i) {
        if(!hanging_node_constraints.is_constrained(i)){
            err_res(i)=system_rhs(i);
           }

        }
    error_residual.u = err_res.l2_norm();
*/

    FullMatrix<double> A;
    A.copy_from(tangent_matrix);
    Vector<double> res(dof_handler.n_dofs());
    double number ;
    number =  A.residual(res,newton_update,system_rhs);// |res| = sys_rhs - A*new_upd
    error_residual.u = number;


}

template <int dim>
void ElasticProblem<dim>::make_constraints(unsigned int &itr){
    std::cout<<" CST "<<std::flush;

    if(itr>=1){
        return;
    }

    hanging_node_constraints.clear();
    DoFTools::make_hanging_node_constraints (dof_handler,
                                             hanging_node_constraints);

    const FEValuesExtractors::Vector disp(0);
    VectorTools::interpolate_boundary_values (dof_handler,
                                              0,
                                              Functions::ZeroFunction<dim>(dim),
                                              hanging_node_constraints
                                              ,fe.component_mask(disp));

    hanging_node_constraints.close ();
}


template <int dim>
Vector<double>
ElasticProblem<dim>::get_total_solution(const Vector<double> &solution_delta) const
{
    Vector<double> current_solution(solution);
    current_solution += solution_delta;
    return current_solution;
}

template <int dim>
void ElasticProblem<dim>::run(AllParameters param){

    using namespace dealii;
/*
    Timer timer;
    timer.start();

    for (unsigned int cy =0; cy<param.fesys.cycles; ++cy)
      {
        std::cout << "Cycle " << cy << ':' << std::endl;
        if (cy == 0)
          {
            //GridGenerator::hyper_cube (triangulation, -1, 1);
            import_mesh(param);
           }
        else
          refine_grid (param);
        std::cout << "   Number of active cells:       "
                  << triangulation.n_active_cells()
                  << std::endl;
        setup_system ();
        std::cout << "   Number of degrees of freedom: "
                  << dof_handler.n_dofs()
                  << std::endl;
        assemble_system ();
        solve (param);
        output_results (cy);
        std::cout<<"time:"<<timer()<<" solution.norm():"<<solution.l2_norm()<<std::endl;
      }
*/

    for (current_time_step=1;current_time_step<=param.fesys.n_time_steps;++current_time_step) {
        if (current_time_step == 1)
          {
            import_mesh(param);
            triangulation.refine_global (param.fesys.gl_ref);
            setup_system ();
          }
//        else
//          refine_grid (param);
        solution_delta = 0.0;
        solve_nonlinear_newton(param,solution_delta);
        solution = solution_delta;
        output_results(current_time_step);
        std::cout<<" solution.norm():"<<solution.l2_norm()<<std::endl;
        std::cout << "   Number of active cells:       "
                          << triangulation.n_active_cells()
                          << std::endl;
        std::cout << "   Number of degrees of freedom: "
                          << dof_handler.n_dofs()
                          << std::endl;
    }
}

template <int dim>
void ElasticProblem<dim>::make_grid()
{
    HyperCubeWithRefinedHole::generate_grid<dim>(triangulation,2,2,5,6);

}

template class thesis::ElasticProblem<2>;
template class thesis::ElasticProblem<3>;


