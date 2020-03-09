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
#include <deal.II/dofs/dof_renumbering.h>

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
      //friend class Others<dim>;
    private:
      void setup_system ();

      std::pair<unsigned int,double> solve_linear_sys (parameters::AllParameters param,dealii::Vector<double> & newton_update);
      void refine_grid (parameters::AllParameters param);
      void output_results (const unsigned int cycle) const;
      void import_mesh(parameters::AllParameters param);

      void solve_nonlinear_newton(parameters::AllParameters param
                                  ,dealii::Vector<double> &solution_delta);
      dealii::Vector<double> get_total_solution(const dealii::Vector<double> &solution_delta)const;

      dealii::Triangulation<dim>   triangulation;
      const dealii::FESystem<dim>        fe;
      dealii::DoFHandler<dim>      dof_handler;

      dealii::ConstraintMatrix     hanging_node_constraints;
      const dealii::QGauss<dim>    quadrature_formula;

      dealii::SparsityPattern      sparsity_pattern;
      dealii::SparseMatrix<double> tangent_matrix;
      dealii::Vector<double>       solution;
      dealii::Vector<double>       system_rhs;

      dealii::Vector<double>       solution_delta;
      unsigned int                 current_time_step;

      static const dealii::SymmetricTensor<4, dim> stress_strain_tensor;

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
              //member variables
              double u;
          };

      Error/*<dim>*/ error_residual, error_residual_0, error_residual_norm;
      void get_error_residual(Error & error_residual,
                              dealii::Vector<double> & newton_update);
      void make_constraints(unsigned int &itr);
      void assemble_system_matrix (Error & error_residual,
                                   dealii::Vector<double> & newton_update,
                                   unsigned int &new_ite,
                                   bool error);
      void assemble_neumann_forces(Error & error_residual,parameters::AllParameters param,
                                   unsigned int &new_ite
                                   ,bool error);


      void print_header();
      void print_footer();
      void make_grid();


    };

}
