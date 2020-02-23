#include "../include/ElasticProblem.h"
#include "../include/others.h"

using namespace dealii;
using namespace thesis;
using namespace parameters;

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

  std::vector<Vector<double> >      rhs_values (n_q_points,Vector<double>(dim+1));
  const FEValuesExtractors::Vector disp(0);
  const Others<dim> right_hand_side;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell_matrix = 0;
      cell_rhs = 0;

      fe_values.reinit(cell);

      right_hand_side.vector_value_list(fe_values.get_quadrature_points(), rhs_values);

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
                const unsigned int
                component_i = fe.system_to_component_index(i).first;
                cell_rhs(i) += fe_values.shape_value(i,q) *
                               rhs_values[q][component_i] *
                               fe_values.JxW(q);
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

template <int dim>
void ElasticProblem<dim>::solve (AllParameters param)
{

  SolverControl           solver_control (param.fesys.steps ,param.fesys.tol);
  SolverCG<>              cg (solver_control);
  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, param.fesys.relax_prm);
  cg.solve (system_matrix, solution, system_rhs,
            preconditioner);
  hanging_node_constraints.distribute (solution);
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
    grid_in.read_ucd(input_file/*, format*/);

}

template <int dim>
void ElasticProblem<dim>::run(AllParameters param){

    using namespace dealii;

    Timer timer;
    timer.start();

    for (unsigned int cy =0; cy<param.fesys.cycles; ++cy)
      {
        std::cout << "Cycle " << cy << ':' << std::endl;
        if (cy == 0)
          {
            GridGenerator::hyper_cube (triangulation, -1, 1);
            //import_mesh(parameter);
            triangulation.refine_global (param.fesys.gl_ref);
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
}

template class thesis::ElasticProblem<2>;


