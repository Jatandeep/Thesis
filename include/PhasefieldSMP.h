#ifndef PHASEFIELDSMP_H
#define PHASEFIELDSMP_H

#include <deal.II/base/timer.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/function.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_point_data.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/data_out_base.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/solver_selector.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>


#include <deal.II/grid/tria.h>
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

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/component_mask.h>

#include "parameter.h"
#include "others.h"

#include <fstream>
#include <iostream>
#include <string>

namespace thesis
{
  
    /*Class to deal with History function responsible for driving d at present time step*//*Thesis_report:Equation 2.42*/
    class PointHistory
    {
      public:
      /*Constructor*/
      PointHistory(){}
       /*Destructor*/
      ~PointHistory() =default;
      /*mutable-> because of need to modify this variable called by a vector of const pointers(CellDataStorage::get_data()const)*/
      double mutable history;
    };

    /*Main class responsible for solving phase field formulation (consisting of coupled PDE) for displacement and phase field.*/
    template <int dim>
    class Phasefield
    {
    public:
      /*Constructor:Takes argument a parameter object "param", contaning all info extracted from input parameter file*/
      Phasefield (const parameters::AllParameters &param);
      /*Destructor*/
      ~Phasefield ();
      /*Data member function containing time loop that dispatches all the work to private functions of this class
      and takes parameter object("param") for passing on info of input parameter file and name of the file("filename") for output files of .vtu format
      and statistics file.*/
      void run (const parameters::AllParameters &param,const std::string filename);

    private:
      /*Forward declaration of a number of objects that are used in parallelizing work using the WorkStream object.
      Declared inside this class as object of these struct(PerTaskData_d/u) will be used to transfer info from
      local to global vectors in function distribute_local_to_global(..). FEValues object are hold in the object
      from ScratchData_d/u which will be used in assembly process. */
      struct PerTaskData_d;
      struct ScratchData_d;
      struct PerTaskData_u;
      struct ScratchData_u;

      /*!Read mesh from external ABAQUS file, perform scaling of the mesh, do global refinement(if specified),
      write .vtk output of input mesh and set boundary ids of read mesh. Takes parameter object "param" as argument
      to do all the above functions.*/
      void import_mesh(const parameters::AllParameters &param);

      /*!Distribute dof's based on a given Finite Element space and allocating memory for the
       * sparse matrix and all used vectors. Set hanging node constraints if local refinement is performed.*/
      void setup_system ();

      /*!Newton-Raphson algorithm looping over all newton iterations for both d and u. Takes parameter object "param"
      for passing information like pre-existing crack strategy, max newton iteration allowed etc. Blockvector
      "solution_delta" contains all the updated solution info for d and u. "delta_t" argument provides updated info
      of delta_t to assembly process of d.(If delta_t is changed after "Time change point"). This function
      returns new iterations for u to check if time adaptivity is required or not in run(..).*/
      unsigned int  solve_nonlinear_newton(const parameters::AllParameters &param
                                  ,dealii::BlockVector<double> &solution_delta
                                  ,const double delta_t);

      /*!Solve the linear system as assembled for d taking parameter object "param" for info like relaxation parameter
      for cg solver. Information of solution updation for a particular newton iteration for d is 
      stored in "newton_update". This function returns a pair linear iterations and linear residual take by
      cg solver.*/
      std::pair<unsigned int,double> solve_sys_d (const parameters::AllParameters &param,
                                                       dealii::BlockVector<double> & newton_update);

      /*!Solve the linear system as assembled for u. Information of solution updation for a particular newton iteration
      for u is stored in "newton_update". This function returns a pair linear iterations and linear residual take by
      SparseDirectUMFPACK which is 0,0 (option provided if any other solver needs to be used).*/
      std::pair<unsigned int,double> solve_sys_u (dealii::BlockVector<double> & newton_update);

      /*!Set hanging node constraints and apply Dirichlet bc for u. Takes number of newton iterations("itr")
      as input as inhomogeneous bc needs to be applied in only first iteration. "load_ratio" decides the 
      per step displacement load to be applied. "param" is used to deduce options like pre-existing crack
      strategy.*/
      void make_constraints_u(unsigned int &itr,const double load_ratio,const parameters::AllParameters &param);

      /*Generate d along slit for M_Id*.Newton Iteration("itr") required as we only need to provide Dirichlet condition
      of d=1 in first iteration. "param" required to identify test case of tension.(Provided if extended for
       shear case in future.)*/
      void make_constraints_d(unsigned int &itr,const parameters::AllParameters &param);

      /*!Assemble the linear system for the d in SMP format. "delta_t" argument provides updated info
      of delta_t to assembly process of d.(If delta_t is changed after "Time change point"). "update" provides
      info to assemble the system at current solution. "param" provides info like material parameters lambd and mu
      for assembly process in assemble_system_d_one_cell(..).*/
      void assemble_system_d (const parameters::AllParameters &param,
                                    dealii::BlockVector<double> & update
                                    ,const double delta_t); 

      /*Wrapper function for d that is executed to do the work in the WorkStream model on one cell for d."param" is passing info
      to assembly process,"cell" contains info for that particular cell. "scratch" and "data" facilitates Workstream functioning.
      For more details about these pls see Step-44 of dealii tutorial.*/                              
      void assemble_system_d_one_cell (const parameters::AllParameters &param,
					                            const typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
                                     	ScratchData_d &scratch,
                                     	PerTaskData_d &data)const;

      /*Wrapper function for d that copies the work done on one cell into the global object that represents it.*/
      void copy_local_to_global_d(const PerTaskData_d &data);


      /*!Assemble the linear system for the u in SMP format."param" provides info like material parameters lambd and mu
      for assembly process in assemble_system_u_one_cell(..)."update" provides info to assemble the system at current solution.*/
      void assemble_system_u(const parameters::AllParameters &param,
			                       dealii::BlockVector<double> & update);

      /*Wrapper function for u that is executed to do the work in the WorkStream model on one cell for d."param" is passing info
      to assembly process,"cell" contains info for that particular cell. "scratch" and "data" facilitates Workstream functioning.
      For more details about these pls see Step-44 of dealii tutorial.*/
      void assemble_system_u_one_cell (const parameters::AllParameters &param,
				                              const typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
                                      ScratchData_u &scratch,
                                      PerTaskData_u &data)const;

      /*Wrapper function for u that copies the work done on one cell into the global object that represents it.*/
      void copy_local_to_global_u(const PerTaskData_u &data);

   
       /*Print header and footer for newton iterations for d and u*/
      void print_header_d();
      void print_footer_d();
      void print_header_u();
      void print_footer_u();

      /*!Write visualization files into .vtu format and generates a solution.pvd file for overall visualization.
      "param" is used for controling output frequency."cycle" is the current timestep number. "filename" is used
      to generate .vtu files with same name as the input parameter file.*/
      void output_results (const parameters::AllParameters &param,unsigned int cycle,const std::string filename) const;

      /*Data memeber variables for implementing the discretization*/
      dealii::Triangulation<dim>   triangulation_m;
      const dealii::FESystem<dim>  fe_m;
      dealii::DoFHandler<dim>      dof_handler_m;

      dealii::ConstraintMatrix     constraints_m;
      const dealii::QGauss<dim>    quadrature_formula_m;
      const dealii::QGauss<dim-1>  face_quadrature_formula_m;

      /*Implementing Block structure to deal with coupled problem of d and u*/
      dealii::BlockSparsityPattern      sparsity_pattern_m;
      dealii::BlockSparseMatrix<double> tangent_matrix_m;
      dealii::BlockVector<double>       solution_m;
      dealii::BlockVector<double>       system_rhs_m;
      /*For time discretization of phasefield, we introduce old solution*/
      dealii::BlockVector<double>       old_solution_m;

      /*To extract d and u part from Block vectors. u varies from 0 to dim-1. d->dim*/
      dealii::FEValuesExtractors::Vector u_extractor;
      dealii::FEValuesExtractors::Scalar d_extractor;

      /*!A struct used to keep track of data needed as convergence criteria.*/
      struct Error
          {
              /*Constructor with data members for calculating error for d and u */
              Error()
              :
              norm(1.0),u(1.0),d(1.0)
              {}
              /*To reset the values after each time step*/
              void reset()
              {
		            norm = 1.0;
                u = 1.0;
		            d = 1.0;
              }
              /*To calculate the normalized version*/
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
      /*Objects with suffix "_residual" represents residual equation error. "_update" suffix represents newton update
      error. "_0" suffix denotes value at first newton iteration. "_norm" suffix denotes the normalized version.*/
      Error error_residual, error_residual_0, error_residual_norm,
	    error_update, error_update_0, error_update_norm;

      /*Calculate error residual from system_rhs for d."error_residual" contains info about residual equation error for d
      which is calculated inside this function.*/
      void get_error_residual_d(Error & error_residual);

      /*Calculate newton error for d. "error_update" contains the error info for d for that particular newton iteration.
      "newton_update" contains updation info for d in that particular newton iteration.*/
      void get_newton_update_error_d(const dealii::BlockVector<double> &newton_update
		      		                    ,Error & error_update);

      /*Calculate error residual from system_rhs for u."error_residual" contains info about residual equation error for u
      which is calculated in this function.*/
      void get_error_residual_u(Error & error_residual);

      /*Calculate newton error for u. "error_update" contains the error info for u for that particular newton iteration.
      "newton_update" contains updation info for u in that particular newton iteration.*/
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

      /*Giving each boundary of the geometry a particular id. If not assigned,default is 0.
      "param" is used to identify "Problem type", whether it is "benchmark" or "lefm_mode_I"*/
      void set_boundary_id(const parameters::AllParameters &param);

      /*dealii class for handing history variables. A class for storing at each cell represented 
      by iterators of type "active_cell_iterator" a vector of "PointHistory" as declared earlier.*/
      dealii::CellDataStorage<typename dealii::Triangulation<dim>::active_cell_iterator
                              ,PointHistory> quadrature_point_history;

      /*Setting up quadrature point history for facilitating history function implementation and 
      initializing data with help of initialize() method of CellDataStorage class.*/
      void setup_qph();

      /*Data member function for calculating the history function for present time step*//*Thesis_report:Equation 2.42*/
      /*"lambda" and "mu" are lame constants, "eps" is strain->all passed for calculating history variable
      which is returned as double value.*/
      double get_history(const double lambda
                	      ,const double mu
    			              ,const dealii::SymmetricTensor<2,dim> &eps)const;

      /*Generating load on boundaries for being written to statistics output file. 
      "param" is used to identify "Problem type", whether it is "benchmark" or "lefm_mode_I"
      "solution" is present time solution used to extract d(for degradation function) and strain
      (for calculating positive and negative stress)*/
      void compute_load(const parameters::AllParameters &param,dealii::BlockVector<double> &solution); 
      
      /*Additional functional to generate a visualization file when load becomes 0 for tension case.
      Takes same argument as compute_load(..). For "tension" case, it returns a load value on top boundary in y-direction which is
      compared to a value of 1. If less than 1, a output visualization file is generated. For other than 
      "tension" cases returned value(1) always generate a false condition(1<1) thereby not affecting the program.*/
      double compute_end_load(const parameters::AllParameters &param,dealii::BlockVector<double> &solution); 
      
      /*Genrating elastic and fracture energy values for statistics file."param" is used to extract 
      material properties and phase field model properties."update" contains the present time solution
      and is used to extract values of d and strain.*/
      void get_energy_v(const parameters::AllParameters &param,
                            dealii::BlockVector<double> & update);
      
      /*For P_I method, this function prescribe d=1 values on selected nodes for initial time step and maintains them throughout simulation*/
      /*Thesis_report:Section4.3.3*//*Takes number of newton iterations(itr) of d as input, as these conditions needs to be applied in only first iteration.
      "param" is needed to check the initial time.*/
      void get_constrained_initial_d(unsigned int itr
                                    ,const parameters::AllParameters &param);
      
      /*For P_I, this function extracts global node ids from cells where crack is to be placed*//*Thesis_report:Section4.3.3*/
      /*"min_cell_dia" contains the element size which is used to  extract the node ids of elements lying in a particular area.
      "param" is used to check various working strategies.*/   
      void extract_initialcrack_d_index(const double min_cell_dia,const parameters::AllParameters &param);
      
      /*Data member variable to store global node ids for the case of P_I.*/
      std::vector<double> global_index_m;

      /*Data members for representing current time */
      double	          current_time_m;

      /*Quantity for display time spend in each section at the end of the program.*/
      mutable dealii::TimerOutput	timer;
      
      /*To print load and energies during the simulation*/
      dealii::TableHandler		statistics;
    };
}

#endif
