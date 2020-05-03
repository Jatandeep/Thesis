#pragma once
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>
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
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include "parameter.h"
#include "others.h"
#include <deal.II/base/timer.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/function.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <fstream>
#include <iostream>

namespace thesis
{
    template <int dim>
    class Phasefield
    {
    public:
      Phasefield (const parameters::AllParameters &param);
      ~Phasefield ();
      void run (const parameters::AllParameters &param);

    private:

      /*!Read mesh from external file*/
      void import_mesh(const parameters::AllParameters &param);

      /*!Distribute dof's based on a given Finite Element space and allocating memory for the
       * sparse matrix and all used vectors.*/
      void setup_system ();

      /*Implement adaptive refinement scheme*/
      void refine_grid (const parameters::AllParameters &param);

      /*!Newton-Raphson algorithm looping over all newton iterations*/
      void solve_nonlinear_newton(const parameters::AllParameters &param
                                  ,dealii::BlockVector<double> &solution_delta);

      /*!Solve the linear system as assembled via assemble_system()*/
      std::pair<unsigned int,double> solve_linear_sys (const parameters::AllParameters &param,
                                                       dealii::BlockVector<double> & newton_update);

      /*!Set hanging node and apply Dirichlet bc.*/
      void make_constraints(unsigned int &itr,const double time);

      /*!Assemble the linear system for the elasticity problem*/
      void assemble_system (const parameters::AllParameters &param,
                            dealii::BlockVector<double> & newton_update);

      /*Assemble External forces(body forces + Neumann forces)*/
      void assemble_rhs(const parameters::AllParameters &param,
			       dealii::BlockVector<double> & newton_update);
     
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
      const dealii::QGauss<dim-1>  face_quadrature_formula_m;

      dealii::BlockSparsityPattern      sparsity_pattern_m;
      dealii::BlockSparseMatrix<double> tangent_matrix_m;
      dealii::BlockVector<double>       solution_m;
      dealii::BlockVector<double>       system_rhs_m;
      dealii::BlockVector<double>       old_solution_m;//for time discretization of phasefield

      double                       current_time_m;

      /*!A struct used to keep track of data needed as convergence criteria.*/
      struct Error
          {
              Error()
              :
              norm(1.0),u(1.0),d(1.0)
              {}

              void reset()
              {
		  norm = 1.0;
                  u = 1.0;
		  d = 1.0;
              }
              void normalize(const Error &err)
              {
                  if (err.norm != 0.0)
                    norm /= err.norm;
                  if (err.u != 0.0)
                    u /= err.u;
                  if(err.d !=0.0)
		    d /= err.d;
              }
              double norm,u,d;
          };

      Error error_residual, error_residual_0, error_residual_norm;

      /*Calculate error residual from system_rhs*/
      void get_error_residual(Error & error_residual);

      dealii::FEValuesExtractors::Vector u_extractor;
      dealii::FEValuesExtractors::Scalar d_extractor;
	
      static const unsigned int n_blocks_m = 2;
      static const unsigned int n_components_m = dim+1;
      static const unsigned int u_dof_start_m = 0;
      static const unsigned int d_dof_start_m = dim;
          
      enum
      {
      u_dof_m = 0,
      d_dof_m = 1
      };
	
      std::vector<dealii::types::global_dof_index> dofs_per_block_m;
      std::vector<dealii::types::global_dof_index> u_element_indices_m;
      std::vector<dealii::types::global_dof_index> d_element_indices_m;

      std::vector<unsigned int> dof_block_identifier_m;
      std::vector<unsigned int> n_comp_per_block{std::vector<unsigned int>(dim,1)};

      void determine_comp_extractor();
      void set_boundary_id();

      double history_m;
      double old_history_m = 0;     
//      double prev_error_norm_m =0;
//      std::shared_ptr<typename InnerPreconditioner<dim>::type/*SparseDirectUMFPACK*/> A_preconditioner;

    };

}
