#pragma once
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/grid/grid_in.h>
#include "parameter.h"
#include "others.h"
#include <deal.II/base/timer.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/function.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/numerics/solution_transfer.h>

#include <fstream>
#include <iostream>


namespace thesis
{
    template <int dim>
    class ElasticProblem
    {
    public:
      ElasticProblem (parameters::AllParameters param);
      virtual ~ElasticProblem ();
      void run (parameters::AllParameters param);

    private:

      /*!Read mesh from external file*/
      void import_mesh(parameters::AllParameters param);

      /*!Distribute dof's based on a given Finite Element space and allocating memory for the
       * sparse matrix and all used vectors.*/
      void setup_system ();

      /*Implement adaptive refinement scheme*/
      void refine_grid (parameters::AllParameters param);

      /*!Newton-Raphson algorithm looping over all newton iterations*/
      void solve_nonlinear_newton(parameters::AllParameters param
                                  ,dealii::Vector<double> &solution_delta);

      /*!Solve the linear system as assembled via assemble_system()*/
      std::pair<unsigned int,double> solve_linear_sys (parameters::AllParameters param,
                                                       dealii::Vector<double> & newton_update);

      /*!Set hanging node and Dirichlet constraints.*/
      void make_constraints(unsigned int &itr);

      /*!Assemble the linear system for the elasticity problem*/
      void assemble_system (parameters::AllParameters param,
                            dealii::Vector<double> & newton_update);

      /*Assemble External forces(body forces + Neumann forces)*/
      void assemble_body_forces(parameters::AllParameters param);

      /*Print header and footer for newton iterations*/
      void print_header();
      void print_footer();

      /*!Write output into files*/
      void output_results (const double cycle) const;

      dealii::Triangulation<dim>   triangulation_m;
      const dealii::FESystem<dim>  fe_m;
      dealii::DoFHandler<dim>      dof_handler_m;

      dealii::ConstraintMatrix     constraints_m;
      const dealii::QGauss<dim>    quadrature_formula_m;

      dealii::SparsityPattern      sparsity_pattern_m;
      dealii::SparseMatrix<double> tangent_matrix_m;
      dealii::Vector<double>       solution_m;
      dealii::Vector<double>       system_rhs_m;

      double                       current_time;

      //dealii::FEValuesExtractors::Vector &disp_extractor(0);

      /*!A struct used to keep track of data needed as convergence criteria.*/
      struct Error
          {
              Error()
              :
              u(1.0)
              {}

              void reset()
              {
                  u = 1.0;
              }
              void normalize(const Error &err)
              {
                  if (err.u != 0.0)
                  {
                      u /= err.u;
                  }
              }
              double u;
          };

      Error error_residual, error_residual_0, error_residual_norm;

      /*Calculate error residual from system_rhs*/
      void get_error_residual(Error & error_residual);

      dealii::SymmetricTensor<2,dim> get_stress(parameters::AllParameters param,dealii::Vector<double> &solution,
                                                const unsigned int q, const unsigned int n_q_points,
                                                dealii::FEValues<dim> &fe_values,
                                                dealii::FEValuesExtractors::Vector disp);
    };

}
