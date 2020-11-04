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
#include <string>

namespace thesis
{
    /*Class to deal with History function responsible for driving d at present time step*/
    class PointHistory
    {
      public:
      PointHistory(){}
      ~PointHistory() =default;
      double mutable history;
    };

    /*Main class responsible for solving phase field formulation*/
    template <int dim>
    class Phasefield
    {
    public:
      /*Constructor*/
      Phasefield (const parameters::AllParameters &param);
      /*Destructor*/
      ~Phasefield ();
      /*Data member function containing time loop*/
      void run (const parameters::AllParameters &param,const std::string filename);

    private:
      /*Variables for dealing with SMP formulation for d and u respectively*/
      struct PerTaskData_d;
      struct ScratchData_d;
      struct PerTaskData_u;
      struct ScratchData_u;

      /*!Read mesh from external ABAQUS file*/
      void import_mesh(const parameters::AllParameters &param);

      /*!Distribute dof's based on a given Finite Element space and allocating memory for the
       * sparse matrix and all used vectors.*/
      void setup_system ();

      /*!Newton-Raphson algorithm looping over all newton iterations for both d and u*/
      unsigned int  solve_nonlinear_newton(const parameters::AllParameters &param
                                  ,dealii::BlockVector<double> &solution_delta
                                  ,const double delta_t);

      /*!Solve the linear system as assembled for d*/
      std::pair<unsigned int,double> solve_sys_d (const parameters::AllParameters &param,
                                                       dealii::BlockVector<double> & newton_update);
      /*!Solve the linear system as assembled for u*/
      std::pair<unsigned int,double> solve_sys_u (dealii::BlockVector<double> & newton_update);

      /*!Set hanging node constraints and apply Dirichlet bc for u.*/
      void make_constraints_u(unsigned int &itr,const double load_ratio,const parameters::AllParameters &param);

      /*Generate d along slit for M_Id*/
      void make_constraints_d(unsigned int &itr,const parameters::AllParameters &param);

      /*!Assemble the linear system for the d in SMP format*/
      void assemble_system_d (const parameters::AllParameters &param,
                                    dealii::BlockVector<double> & update
                                    ,const double delta_t);  
      void assemble_system_d_one_cell (const parameters::AllParameters &param,
					                            const typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
                                     	ScratchData_d &scratch,
                                     	PerTaskData_d &data)const;
      void copy_local_to_global_d(const PerTaskData_d &data);


      /*!Assemble the linear system for the u in SMP format*/
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
      void output_results (const parameters::AllParameters &param,unsigned int cycle,const std::string filename) const;

      /*Data memeber variables for implementing the discretization*/
      dealii::Triangulation<dim>   triangulation_m;
      const dealii::FESystem<dim>  fe_m;
      dealii::DoFHandler<dim>      dof_handler_m;

      dealii::ConstraintMatrix     constraints_m;
      const dealii::QGauss<dim>    quadrature_formula_m;
      const dealii::QGauss<dim-1>  face_quadrature_formula_m;

      /*Implementing Block structure to delai with coupled problem of d and u*/
      dealii::BlockSparsityPattern      sparsity_pattern_m;
      dealii::BlockSparseMatrix<double> tangent_matrix_m;
      dealii::BlockVector<double>       solution_m;
      dealii::BlockVector<double>       system_rhs_m;
      /*For time discretization of phasefield, we introduce old solution*/
      dealii::BlockVector<double>       old_solution_m;

      /*To extract d and u part from Block vectors*/
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

      /*Calculate error residual from system_rhs for d*/
      void get_error_residual_d(Error & error_residual);
      /*Calculate newton error for d*/
      void get_newton_update_error_d(const dealii::BlockVector<double> &newton_update
		      		                    ,Error & error_update);

      /*Calculate error residual from system_rhs for d*/
      void get_error_residual_u(Error & error_residual);
      /*Calculate newton error for u*/
      void get_newton_update_error_u(const dealii::BlockVector<double> &newton_update
		      		                    ,Error & error_update);

      /*Defining basic structure for dealing with block vectors*/      
      static const unsigned int n_blocks_m = 2;
      static const unsigned int n_components_m = dim+1;
      static const unsigned int u_dof_start_m = 0;
      static const unsigned int d_dof_start_m = dim;

      /*Assigning a dedicated number for recognizing d and u blocks*/    
      enum
      {
      u_dof_m = 0,
      d_dof_m = 1
      };
	
      /*Data member variables for total dofs per block*/
      std::vector<dealii::types::global_dof_index> dofs_per_block_m;
      /*Data member variables for identifying d and u contributions per cell*/
      std::vector<unsigned int> dof_block_identifier_m;

      /*Data member function for identifying the block for d and u respectively*/
      void determine_comp_extractor();
      /*Giving each boundary of the geometry a particular id. If not assigned,default is 0*/
      void set_boundary_id(const parameters::AllParameters &param);

      /*dealii function for handing history variables*/
      dealii::CellDataStorage<typename dealii::Triangulation<dim>::active_cell_iterator
                              ,PointHistory> quadrature_point_history;
      /*Setting up quadrature point history for facilitating history function implementation*/
      void setup_qph();
      /*Data member function for calculating the history function for present time step*/
      double get_history(const double lambda
                	      ,const double mu
    			              ,const dealii::SymmetricTensor<2,dim> &eps)const;

      /*Generating load on boundaries for being written to statistics file*/
      void compute_load(const parameters::AllParameters &param,dealii::BlockVector<double> &solution); 
      /*Additional functional to generate a visualization file when load becomes 0 for tension case*/
      double compute_end_load(const parameters::AllParameters &param,dealii::BlockVector<double> &solution); 
      /*Genrating elastic and fracture energy values for statistics file*/
      void get_energy_v(const parameters::AllParameters &param,
                            dealii::BlockVector<double> & update);
      /*For P_I method, this function prescribe d=1 values on selected nodes for initial time step and maintains them throughout simulation*/
      void get_constrained_initial_d(unsigned int itr
                                    ,const parameters::AllParameters &param);
      /*For P_I, this function extracts global node ids from cells where crack is to be placed*/   
      void extract_initialcrack_d_index(const double min_cell_dia,const parameters::AllParameters &param);
      /*Data member varaible to store global node ids*/
      std::vector<double> global_index_m;

      double	          current_time_m;
      mutable dealii::TimerOutput	timer;
      /*To print load and energies during the simulation*/
      dealii::TableHandler		statistics;
    };
}
