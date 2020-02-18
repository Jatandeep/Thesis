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

#include <fstream>
#include <iostream>



    using namespace dealii;

    template <int dim>
    class ElasticProblem
    {
    public:
      ElasticProblem (const std::string &filename);
      virtual ~ElasticProblem ();
      void run ();

    private:
      void setup_system ();
      void assemble_system ();
      void solve ();
      void refine_grid ();
      void output_results (const unsigned int cycle) const;

      Triangulation<dim>   triangulation;
      AllParameters parameter;
      const FESystem<dim>        fe;

      DoFHandler<dim>      dof_handler;

      ConstraintMatrix     hanging_node_constraints;
      const QGauss<dim>    quadrature_formula;
      SparsityPattern      sparsity_pattern;
      SparseMatrix<double> system_matrix;
      Vector<double>       solution;
      Vector<double>       system_rhs;

      void import_mesh();

      static const SymmetricTensor<4, dim> stress_strain_tensor;

    };













