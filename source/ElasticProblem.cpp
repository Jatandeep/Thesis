#include "../include/ElasticProblem.h"
#include "../include/others.h"

using namespace dealii;
using namespace thesis;
/*
template <int dim>
ElasticProblem<dim>::ElasticProblem(const std::string &filename):
    parameter(filename),fe(FE_Q<dim>(parameter.fe_degree),dim), dof_handler(triangulation),quadrature_formula(parameter.fe_degree+1)
{}
*/
template <int dim>
SymmetricTensor<4, dim> get_stress_strain_tensor(const double lambda,
                                                 const double mu)
{
  SymmetricTensor<4, dim> tmp;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        for (unsigned int l = 0; l < dim; ++l)
          tmp[i][j][k][l] = (((i == j) && (k == l)) ? lambda:0)
                            + (((i == k) && (j == l)) ? mu:0)
                            + (((i == l) && (j == k)) ? mu:0);
  return tmp;
}
template <int dim>
inline SymmetricTensor<2, dim> get_strain(const FEValues<dim> &fe_values,
                                          const unsigned int   shape_func,
                                          const unsigned int   q_point)
{
  SymmetricTensor<2, dim> temp;

  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i ; j < dim; ++j)
      temp[i][j] =
        (fe_values.shape_grad_component(shape_func, q_point, i)[j] +
         fe_values.shape_grad_component(shape_func, q_point, j)[i]) /
        2;
  return temp;
}

template <int dim>
const SymmetricTensor<4, dim> ElasticProblem<dim>::stress_strain_tensor =
  get_stress_strain_tensor<dim>(1, 1);

template <int dim>
ElasticProblem<dim>::ElasticProblem(AllParameters parameter):
    fe(FE_Q<dim>(parameter.fe_degree),dim), dof_handler(triangulation),quadrature_formula(parameter.fe_degree+1)
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
  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
                                           hanging_node_constraints);
  hanging_node_constraints.close ();
  DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler,
                                  dsp,
                                  hanging_node_constraints,
                                  /*keep_constrained_dofs = */ true);
  sparsity_pattern.copy_from (dsp);
  system_matrix.reinit (sparsity_pattern);
  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
}

template <int dim>
void ElasticProblem<dim>::assemble_system ()
{

  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values   | update_gradients |
                           update_quadrature_points | update_JxW_values);
  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();
  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  std::vector<double>     lambda_values (n_q_points);
  std::vector<double>     mu_values (n_q_points);
  Functions::ConstantFunction<dim> lambda(1.), mu(1.);
  std::vector<Tensor<1, dim> > rhs_values (n_q_points);
//  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
//                                                 endc = dof_handler.end();

  const FEValuesExtractors::Vector disp(0);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell_matrix = 0;
      cell_rhs = 0;

      fe_values.reinit(cell);

//Step-8 Definition
/*
      lambda.value_list (fe_values.get_quadrature_points(), lambda_values);
      mu.value_list     (fe_values.get_quadrature_points(), mu_values);
      right_hand_side (fe_values.get_quadrature_points(), rhs_values);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          const unsigned int
          component_i = fe.system_to_component_index(i).first;
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            {
              const unsigned int
              component_j = fe.system_to_component_index(j).first;
              for (unsigned int q_point=0; q_point<n_q_points;
                   ++q_point)
                {
                  cell_matrix(i,j)
                  +=
                    (
                      (fe_values.shape_grad(i,q_point)[component_i] *
                       fe_values.shape_grad(j,q_point)[component_j] *
                       lambda_values[q_point])
                      +
                      (fe_values.shape_grad(i,q_point)[component_j] *
                       fe_values.shape_grad(j,q_point)[component_i] *
                       mu_values[q_point])
                      +
                      ((component_i == component_j) ?
                       (fe_values.shape_grad(i,q_point) *
                        fe_values.shape_grad(j,q_point) *
                        mu_values[q_point])  :
                       0)
                    )
                    *
                    fe_values.JxW(q_point);
                }
            }
        }
*/

        Others<dim>::right_hand_side.value_list(fe_values.get_quadrature_points(), rhs_values);
/*
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
              {
                const SymmetricTensor<2, dim>
                  eps_phi_i = get_strain(fe_values, i, q_point),
                  eps_phi_j = get_strain(fe_values, j, q_point);
                cell_matrix(i, j) += (eps_phi_i *
                                      stress_strain_tensor *
                                      eps_phi_j
                                      ) *
                                     fe_values.JxW(q_point);
              }

        for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
             const unsigned int
                  component_i = fe.system_to_component_index(i).first;
                  for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
                    cell_rhs(i) += fe_values.shape_value(i,q_point) *
                                   rhs_values[q_point][component_i] *
                                   fe_values.JxW(q_point);
             }

        cell->get_dof_indices (local_dof_indices);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
             for (unsigned int j=0; j<dofs_per_cell; ++j)
                          system_matrix.add (local_dof_indices[i],
                                             local_dof_indices[j],
                                             cell_matrix(i,j));
                        system_rhs(local_dof_indices[i]) += cell_rhs(i);
            }
       }
*/

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
                cell_rhs(i) += phi_i_u *
                               rhs_values[q]*
                               fe_values.JxW(q);
            }




       hanging_node_constraints.condense (system_matrix);
       hanging_node_constraints.condense (system_rhs);
       std::map<types::global_dof_index,double> boundary_values;
       VectorTools::interpolate_boundary_values (dof_handler,
                                                          0,
                                                          Functions::ZeroFunction<dim>(dim),
                                                          boundary_values);
       MatrixTools::apply_boundary_values (boundary_values,
                                                    system_matrix,
                                                    solution,
                                                    system_rhs);
}
}

template <int dim>
void ElasticProblem<dim>::solve (AllParameters parameter)
{

  SolverControl           solver_control (parameter.steps ,parameter.tol);
  SolverCG<>              cg (solver_control);
  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, parameter.relax_prm);
  cg.solve (system_matrix, solution, system_rhs,
            preconditioner);
  hanging_node_constraints.distribute (solution);
}

template <int dim>
void ElasticProblem<dim>::refine_grid (AllParameters parameter)
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
  KellyErrorEstimator<dim>::estimate (dof_handler,
                                      QGauss<dim-1>(parameter.quad_order),
                                      typename FunctionMap<dim>::type(),
                                      solution,
                                      estimated_error_per_cell);
  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                   estimated_error_per_cell,
                                                   parameter.act_ref, parameter.act_cors);
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
void ElasticProblem<dim>::import_mesh(AllParameters parameter){

    std::string grid_name;
    grid_name += parameter.meshfile;

    GridIn<dim> grid_in;
    grid_in.attach_triangulation(triangulation);
    std::ifstream input_file(grid_name.c_str());
    Assert(dim==2, ExcInternalError());
    grid_in.read_ucd(input_file/*, format*/);

}

template <int dim>
void ElasticProblem<dim>::run(AllParameters parameter){

    using namespace dealii;

    Timer timer;
    timer.start();

    for (unsigned int cy =0; cy<parameter.cycles; ++cy)
      {
        std::cout << "Cycle " << cy << ':' << std::endl;
        if (cy == 0)
          {
            //GridGenerator::hyper_cube (triangulation, -1, 1);
            import_mesh(parameter);
            triangulation.refine_global (parameter.gl_ref);
          }
        else
          refine_grid (parameter);
        std::cout << "   Number of active cells:       "
                  << triangulation.n_active_cells()
                  << std::endl;
        setup_system ();
        std::cout << "   Number of degrees of freedom: "
                  << dof_handler.n_dofs()
                  << std::endl;
        assemble_system ();
        solve (parameter);
        output_results (cy);
        std::cout<<"time:"<<timer()<<" solution.norm():"<<solution.l2_norm()<<std::endl;
      }
}

template class thesis::ElasticProblem<2>;


