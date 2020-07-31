#pragma once

#include <deal.II/base/timer.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/function.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_point_data.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/table_handler.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/packaged_operation.h>
#include <deal.II/lac/precondition_selector.h>
#include <deal.II/lac/solver_selector.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/solution_transfer.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/mapping_q_eulerian.h>

#include "parameter.h"
#include "others.h"

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
      unsigned int solve_nonlinear_newton(const parameters::AllParameters &param
                                  ,dealii::BlockVector<double> &solution_delta
                                  ,double delta_t);

      /*!Solve the linear system as assembled via assemble_system()*/
      std::pair<unsigned int,double> solve_sys_d (const parameters::AllParameters &param,
                                                       dealii::BlockVector<double> & newton_update);
                                                     /*dealii::Vector<double> & newton_update);*/

      /*!Solve the linear system as assembled via assemble_system()*/
      std::pair<unsigned int,double> solve_sys_u (const parameters::AllParameters &param,
                                                       dealii::BlockVector<double> & newton_update);

      /*!Set hanging node and apply Dirichlet bc.*/
      void make_constraints(unsigned int &itr,const double load_ratio,const double u_total);

      /*!Assemble the linear system for the elasticity problem*/
      void assemble_system_d (const parameters::AllParameters &param,
                            dealii::BlockVector<double> & update
                          /*,dealii::Vector<double> & update*/
                            ,double delta_t);
      
      /*!Assemble the linear system for the elasticity problem*/
      void assemble_system_u (const parameters::AllParameters &param,
                            dealii::BlockVector<double> & update);

      /*Print header and footer for newton iterations*/
      void print_header_d();
      void print_footer_d();
      void print_header_u();
      void print_footer_u();

      /*!Write output into files*/
      void output_results (const parameters::AllParameters &param,unsigned int cycle) const;

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

      dealii::FEValuesExtractors::Vector u_extractor;
      dealii::FEValuesExtractors::Scalar d_extractor;
	

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

      Error error_residual, error_residual_0, error_residual_norm,
	    error_update, error_update_0, error_update_norm;

      /*Calculate error residual from system_rhs*/
      void get_error_residual_d(Error & error_residual);
      void get_newton_update_error_d(const dealii::BlockVector<double> &newton_update
                                   /*const dealii::Vector<double> &newton_update*/
		      		                      ,Error & error_update);
      
      void get_error_residual_u(Error & error_residual);
      void get_newton_update_error_u(const dealii::BlockVector<double> &newton_update
		      		                    ,Error & error_update);
      
      
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

      /*History class, variables and functions*/

      const dealii::FE_DGQ<dim>  history_fe_m;
      dealii::DoFHandler<dim>    history_dof_handler_m;
      dealii::ConstraintMatrix   history_constraints_m;
      
      struct PointHistory
      {
	      double history;
      };

      void setup_quadrature_point_history(); 
      std::vector<PointHistory> quadrature_point_history;

      std::vector< std::vector< dealii::/*Block*/Vector<double> > >
             history_field {std::vector< std::vector< dealii::/*Block*/Vector<double> > >
				(/*dim*/1, std::vector< dealii::/*Block*/Vector<double> >(/*dim*/1)) },
             local_history_values_at_qpoints {std::vector< std::vector< dealii::/*Block*/Vector<double> > >
				(/*dim*/1, std::vector< dealii::/*Block*/Vector<double> >(/*dim*/1)) },
             local_history_fe_values {std::vector< std::vector< dealii::/*Block*/Vector<double> > >
				(/*dim*/1, std::vector< dealii::/*Block*/Vector<double> >(/*dim*/1)) };
     
      void history_quadrature_to_global();
      void history_global_to_quadrature();
    
      double get_history(const double lambda
                        ,const double mu
                        ,const dealii::SymmetricTensor<2,dim> &eps)const;

      void compute_load(const double lambda,const double mu,dealii::BlockVector<double> &solution);    
      void project_back_phase_field();

      double get_energy(const parameters::AllParameters &param,
                            dealii::BlockVector<double> & update);
      
      void compute_lefm_errors(const parameters::AllParameters &param);
      double		                current_time_m;
      mutable dealii::TimerOutput	timer;
      dealii::TableHandler              statistics;
    };
}