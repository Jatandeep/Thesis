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
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/index_set.h>

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
#include <deal.II/dofs/block_info.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/mapping_q_eulerian.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/component_mask.h>

#include "parameter.h"
#include "others.h"

#include <fstream>
#include <iostream>


namespace thesis
{
    
    class PointHistory
    {
      public:
      PointHistory(){}
      ~PointHistory() =default;
      double mutable history;
    };

    template <int dim>
    class Phasefield
    {
    public:
      Phasefield (const parameters::AllParameters &param);
      ~Phasefield ();
      void run (const parameters::AllParameters &param);

    private:

      struct PerTaskData_d;
      struct ScratchData_d;

      struct PerTaskData_u;
      struct ScratchData_u;

      /*!Read mesh from external file*/
      void import_mesh(const parameters::AllParameters &param);

      /*!Distribute dof's based on a given Finite Element space and allocating memory for the
       * sparse matrix and all used vectors.*/
      void setup_system ();

      /*Implement adaptive refinement scheme*/
      void refine_grid (const parameters::AllParameters &param);

      /*!Newton-Raphson algorithm looping over all newton iterations*/
      unsigned int  solve_nonlinear_newton(const parameters::AllParameters &param
                                  ,dealii::BlockVector<double> &solution_delta
                                  ,const double delta_t);

      /*!Solve the linear system as assembled via assemble_system()*/
      std::pair<unsigned int,double> solve_sys_d (const parameters::AllParameters &param,
                                                       dealii::BlockVector<double> & newton_update);
      /*!Solve the linear system as assembled via assemble_system()*/
      std::pair<unsigned int,double> solve_sys_u (const parameters::AllParameters &param,
                                                       dealii::BlockVector<double> & newton_update);

      /*!Set hanging node and apply Dirichlet bc.*/
      void make_constraints_u(unsigned int &itr,const double load_ratio,const parameters::AllParameters &param);

      /*!Set hanging node and apply Dirichlet bc.*/
      void make_constraints_d(unsigned int &itr,const parameters::AllParameters &param);

      /*!Assemble the linear system for the elasticity problem*/
      void assemble_system_d (const parameters::AllParameters &param,
                                    dealii::BlockVector<double> & update
                                    ,const double delta_t);
  
      void assemble_system_d_one_cell (const parameters::AllParameters &param,
					                            const typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
                                     	ScratchData_d &scratch,
                                     	PerTaskData_d &data)const;

      void copy_local_to_global_d(const PerTaskData_d &data);


      /*Assemble External forces(body forces + Neumann forces)*/
      void assemble_system_u(const parameters::AllParameters &param,
			                       dealii::BlockVector<double> & update);
      void assemble_system_u_one_cell (const parameters::AllParameters &param,
				                              const typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
                                      ScratchData_u &scratch,
                                      PerTaskData_u &data)const;

      void copy_local_to_global_u(const PerTaskData_u &data);

   
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
		      		                    ,Error & error_update);
      
      void get_error_residual_u(Error & error_residual);
      void get_newton_update_error_u(const dealii::BlockVector<double> &newton_update
		      		                    ,Error & error_update);
      
      void get_error_residual(Error & error_residual);
      
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
      /*
      const dealii::FE_DGQ<dim>  history_fe_adp_m;
      dealii::DoFHandler<dim>    history_dof_handler_adp_m;
      dealii::ConstraintMatrix   history_constraints_adp_m;
      
      struct PointHistory_adp
      {
	      double history_adp;
      };
      
      void setup_quadrature_point_history_adp(); 
      std::vector<PointHistory_adp> quadrature_point_history_adp;
           
      dealii::Vector<double> history_field_adp,local_history_values_at_qpoints_adp,local_history_fe_values_adp;

      void history_quadrature_to_global_adp();
      void history_global_to_quadrature_adp();
      */

      dealii::CellDataStorage<typename dealii::Triangulation<dim>::cell_iterator
                              ,PointHistory> quadrature_point_history;
      void setup_qph();

      double get_history(const double lambda
                	      ,const double mu
    			              ,const dealii::SymmetricTensor<2,dim> &eps)const;

      void compute_load(const parameters::AllParameters &param,dealii::BlockVector<double> &solution);    

      std::pair<double,double> get_energy_p(const parameters::AllParameters &param,
                            dealii::BlockVector<double> & update);
      void get_energy_v(const parameters::AllParameters &param,
                            dealii::BlockVector<double> & update);
      
      void compute_lefm_errors(const parameters::AllParameters &param);
      double get_critical_stress(const parameters::AllParameters &param);

      void get_constrained_initial_d(unsigned int itr
                                    ,const parameters::AllParameters &param);
      void extract_initialcrack_d_index(const double min_cell_dia,const parameters::AllParameters &param);
      std::vector<double> global_index_m;

      double		                current_time_m;
      mutable dealii::TimerOutput	timer;
      dealii::TableHandler		statistics;
    };
}
