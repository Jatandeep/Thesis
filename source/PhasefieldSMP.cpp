
#include "../include/PhasefieldSMP.h"
#include "../include/others.h"
#include "../include/constants.h"
#include "../include/constitutive.h"
#include "../include/utilities.h"

using namespace dealii;
using namespace thesis;
using namespace parameters;

template <int dim>
Phasefield<dim>::Phasefield(const AllParameters &param)
  : fe_m(FESystem<dim>(FE_Q<dim>(param.fesys.fe_degree),dim),1,
	FE_Q<dim>(param.fesys.fe_degree),1)
  , quadrature_formula_m(param.fesys.quad_order)
  , face_quadrature_formula_m(param.fesys.quad_order)
  , u_extractor(u_dof_start_m)
  , d_extractor(d_dof_start_m)
  , dofs_per_block_m(n_blocks_m)
  , history_fe_m(1)						//history	
  , timer(std::cout,TimerOutput::summary,TimerOutput::wall_times)
				
{
    import_mesh(param);
    dof_handler_m.initialize(triangulation_m, fe_m);
    history_dof_handler_m.initialize(triangulation_m,history_fe_m); //history
    setup_system();			
    setup_quadrature_point_history();				//history
    determine_comp_extractor();
}



template <int dim>
Phasefield<dim>::~Phasefield ()
{
  dof_handler_m.clear ();
  history_dof_handler_m.clear();			//history
  statistics.set_auto_fill_mode(true);
}


template <int dim>
void Phasefield<dim>::setup_system ()
{
  timer.enter_subsection("Setup system");

  dof_handler_m.distribute_dofs (fe_m);
  
  DoFRenumbering::block_wise(dof_handler_m); 
  DoFTools::count_dofs_per_block (dof_handler_m, dofs_per_block_m);

  std::cout<<std::endl;
  std::cout << "   Number of active cells:       "
                            << triangulation_m.n_active_cells()
                            << std::endl;
  std::cout << "   Number of degrees of freedom: "
                            << dof_handler_m.n_dofs()
                            << std::endl;
/*
  constraints_m.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler_m,		//For Adaptive Ref
                                           constraints_m);
  constraints_m.close ();
*/
 
  tangent_matrix_m.clear();

  {
  BlockDynamicSparsityPattern dsp(dofs_per_block_m,dofs_per_block_m);
  dsp.collect_sizes();

  DoFTools::make_sparsity_pattern(dof_handler_m
                                  ,dsp
                                  ,constraints_m
                                  ,false);

  sparsity_pattern_m.copy_from (dsp);
  }

  tangent_matrix_m.reinit (sparsity_pattern_m);

  system_rhs_m.reinit (dofs_per_block_m);
  system_rhs_m.collect_sizes();

  solution_m.reinit (dofs_per_block_m);
  solution_m.collect_sizes();

//  std::cout<<"setup solution_d.size: "<<solution_m.block(d_dof_m).size()<<std::endl;

  old_solution_m.reinit (dofs_per_block_m);
  old_solution_m.collect_sizes();

  timer.leave_subsection();
}


template <int dim>
void Phasefield<dim>::determine_comp_extractor()
{
  u_element_indices_m.clear();
  d_element_indices_m.clear();
	
  const unsigned int   dofs_per_cell = fe_m.dofs_per_cell;

   for (unsigned int k = 0; k < dofs_per_cell; ++k)
      {
        const unsigned int k_group = fe_m.system_to_base_index(k).first.first;//Return for shape function index(k) the base element it belongs to.
        if (k_group == u_dof_m)
          u_element_indices_m.push_back(k);
        else if (k_group == d_dof_m)
          d_element_indices_m.push_back(k);
        else
          {
            Assert(k_group <= d_dof_m, ExcInternalError());
          }
      }
	dof_block_identifier_m.clear();
	for(unsigned int k=0;k<dofs_per_cell;++k)
		dof_block_identifier_m.push_back(fe_m.system_to_base_index(k).first.first);
}

template <int dim>
void Phasefield<dim>::set_boundary_id()
{
  double tol_machine = 1e-10;
  typename Triangulation<dim>::active_cell_iterator cell =
      triangulation_m.begin_active(), endc = triangulation_m.end();

  for (;cell != endc; ++cell)
	for (unsigned int f=0;f < GeometryInfo<dim>::faces_per_cell;++f)
          {
            const Point<dim> face_center = cell->face(f)->center();
            if (cell->face(f)->at_boundary())
              {
		//////////////////////For (0,1)x(0,1)
		//left boundary
                if ((face_center[0] < 0.0+tol_machine) && (face_center[0] > 0.0-tol_machine)
                   )
                  cell->face(f)->set_boundary_id(0);
		//right boundary
                else if ((face_center[0] < 1.0+tol_machine) && (face_center[0] > 1.0-tol_machine)
                        )
                  cell->face(f)->set_boundary_id(1);
		// bottom boundary
                else if ((face_center[1] < 0.0+tol_machine) && (face_center[1] > 0.0-tol_machine)
                        )
                  cell->face(f)->set_boundary_id(2);
              	// top boundary
                else if ((face_center[1] < 1.0+tol_machine) && (face_center[1] > 1.0-tol_machine)
                        )
                  cell->face(f)->set_boundary_id(3);
		// lower part of slit
               /* else if ((face_center[0] < 1.0) && (face_center[0] > 0.5) && (face_center[1])
		 	)
                  cell->face(f)->set_boundary_id(4);
		*/
/*
		/////////////////////For (-1,1)x(-1,1) // For Reference solution
		//left boundary
                if ((face_center[0] < -1.0+tol_machine) && (face_center[0] > -1.0-tol_machine)
                   )
                  cell->face(f)->set_boundary_id(0);
		//right boundary
                else if ((face_center[0] < 1.0+tol_machine) && (face_center[0] > 1.0-tol_machine)
                        )
                  cell->face(f)->set_boundary_id(1);
		// bottom boundary
                else if ((face_center[1] < -1.0+tol_machine) && (face_center[1] > -1.0-tol_machine)
                        )
                  cell->face(f)->set_boundary_id(2);
              	// top boundary
                else if ((face_center[1] < 1.0+tol_machine) && (face_center[1] > 1.0-tol_machine)
                        )
                  cell->face(f)->set_boundary_id(3);
*/
	      }
          }

}

template <int dim>
void Phasefield<dim>::setup_quadrature_point_history()
{
//    std::cout << "    Setting up quadrature point data..." << std::endl;
    
    history_dof_handler_m.distribute_dofs(history_fe_m);	//history
/*    
    history_constraints_m.clear ();
    DoFTools::make_hanging_node_constraints (history_dof_handler_m,
                                           history_constraints_m);
    history_constraints_m.close ();
*/
    for (unsigned int i=0; i<1; i++)
      for (unsigned int j=0; j<1; j++)
      {								
	history_field[i][j].reinit(history_dof_handler_m.n_dofs());
	local_history_values_at_qpoints[i][j].reinit(quadrature_formula_m.size());
	local_history_fe_values[i][j].reinit(history_fe_m.dofs_per_cell);
        //std::cout<<"setup_system_history"<<std::endl;
      }

    triangulation_m.clear_user_data();
    {
      std::vector<PointHistory> tmp;
      quadrature_point_history.swap(tmp);
    }
    quadrature_point_history.resize(triangulation_m.n_active_cells() * quadrature_formula_m.size());

    unsigned int history_index = 0;
//    for (auto &cell : triangulation_m.active_cell_iterators()){
    typename Triangulation<dim>::active_cell_iterator cell =
      triangulation_m.begin_active(), endc = triangulation_m.end();

    for (;cell != endc; ++cell)
    {
      cell->set_user_pointer(&quadrature_point_history[history_index]);
      history_index += quadrature_formula_m.size();
    }
    Assert(history_index == quadrature_point_history.size(),ExcInternalError());

}

template <int dim>
double Phasefield<dim>::get_energy(const AllParameters &param,BlockVector<double> & update)
{

  FEValues<dim> fe_values (fe_m, quadrature_formula_m,
                           update_values   | update_gradients |
                           update_quadrature_points | update_JxW_values);

  const unsigned int   n_q_points    = quadrature_formula_m.size();

  double Energy_density=0,Dissipation_energy=0,Total_energy=0;
    
  for (const auto &cell : dof_handler_m.active_cell_iterators())
    {
      fe_values.reinit(cell);

      std::vector<SymmetricTensor<2,dim>> epsilon_vals(n_q_points);
      fe_values[u_extractor].get_function_symmetric_gradients(update,epsilon_vals);
 
      std::vector<double> d_vals(n_q_points);
      fe_values[d_extractor].get_function_values(update,d_vals);
      
      std::vector<Tensor<1,dim>> grad_d(n_q_points);
      fe_values[d_extractor].get_function_gradients(update,grad_d);

      for (unsigned int q = 0; q < n_q_points; ++q){
	
		Energy_density += ( get_deg_func(d_vals[q]) + param.pf.k )*get_energy_density_plus(param.materialmodel.lambda,param.materialmodel.mu,epsilon_vals[q]) 
				+ get_energy_density_minus(param.materialmodel.lambda,param.materialmodel.mu,epsilon_vals[q]);

      		Dissipation_energy += param.pf.g_c*( (0.5/param.pf.l)*d_vals[q]*d_vals[q] + 0.5*param.pf.l* (scalar_product(grad_d[q],grad_d[q])) );
      }
    }
  Total_energy = Energy_density + Dissipation_energy;
  
//  std::cout<<"Energy density:"<<Energy_density<<std::endl;
//  std::cout<<"Dissipation Energy:"<<Dissipation_energy<<std::endl;

return Total_energy;
}


template <int dim>
struct Phasefield<dim>::PerTaskData_K
{
    FullMatrix<double>        cell_matrix;
    std::vector<types::global_dof_index> local_dof_indices;
    PerTaskData_K(const FiniteElement<dim> &fe)
      :
      cell_matrix(fe.dofs_per_cell, fe.dofs_per_cell),
      local_dof_indices(fe.dofs_per_cell)
    {}

    void reset()
    {
      cell_matrix = 0.0;
    }
};

template <int dim>
struct Phasefield<dim>::ScratchData_K
{
    FEValues<dim> fe_values;

    std::vector<SymmetricTensor<2,dim>> epsilon_vals;
    std::vector<SymmetricTensor<2,dim>> old_epsilon_vals;
    std::vector<double> d_vals;
    SymmetricTensor<2,dim> tmp1,tmp2;
    SymmetricTensor<4,dim> BigC_plus;
    SymmetricTensor<4,dim> BigC_minus;
    SymmetricTensor<2,dim> sigma_plus;

    const BlockVector<double> & update;
    const BlockVector<double> & old_solution_m;
    std::vector<double> shape_d;
    std::vector<Tensor<2, dim>> grad_shape_u;
    std::vector<Tensor<1, dim>> grad_shape_d;
    std::vector<SymmetricTensor<2, dim>> sym_grad_shape_u;

    ScratchData_K( const FiniteElement<dim> &fe_cell
                  ,const QGauss<dim> &qf_cell
                  ,const UpdateFlags uf_cell
	          ,const BlockVector<double> &upd
		  ,const BlockVector<double> &old_sol_m)
     	    :fe_values(fe_cell, qf_cell, uf_cell)
		,epsilon_vals(qf_cell.size())
		,old_epsilon_vals(qf_cell.size())
		,d_vals(qf_cell.size())
		,update(upd)
		,old_solution_m(old_sol_m)
		,shape_d(fe_cell.dofs_per_cell)
		,grad_shape_u(fe_cell.dofs_per_cell)
		,grad_shape_d(fe_cell.dofs_per_cell)
		,sym_grad_shape_u(fe_cell.dofs_per_cell)
		{}

    ScratchData_K(const ScratchData_K &Scratch)
      		:fe_values(Scratch.fe_values.get_fe()
                ,Scratch.fe_values.get_quadrature()
                ,Scratch.fe_values.get_update_flags())
      		,epsilon_vals(Scratch.epsilon_vals)
		,old_epsilon_vals(Scratch.old_epsilon_vals)
      		,d_vals(Scratch.d_vals)
		,tmp1(Scratch.tmp1)
		,tmp2(Scratch.tmp2)
		,BigC_plus(Scratch.BigC_plus)
		,BigC_minus(Scratch.BigC_minus)
		,sigma_plus(Scratch.sigma_plus)
		,update(Scratch.update)
		,old_solution_m(Scratch.old_solution_m)
		,shape_d(Scratch.shape_d)
		,grad_shape_u(Scratch.grad_shape_u)
       		,grad_shape_d(Scratch.grad_shape_d)
		,sym_grad_shape_u(Scratch.sym_grad_shape_u)
		{}

    void reset()
    {
    const unsigned int n_q_points = epsilon_vals.size();
    for(unsigned int q=0;q<n_q_points;++q)
    {
	epsilon_vals[q] = 0.0;
    	old_epsilon_vals[q] = 0.0;
    	d_vals[q] = 0.0;
    }

    }
};

template <int dim>
void Phasefield<dim>::assemble_system_tangent (const AllParameters &param,BlockVector<double> & update)
{
  timer.enter_subsection("Assemble tangent matrix");

  tangent_matrix_m = 0.0;

  const UpdateFlags uf_cell(update_values    |  update_gradients |
                            update_quadrature_points |
			      update_JxW_values);

  PerTaskData_K per_task_data(fe_m);
  ScratchData_K scratch_data(fe_m, quadrature_formula_m, uf_cell,update,old_solution_m);

  WorkStream::run (dof_handler_m.begin_active(),
                   dof_handler_m.end(),
                   std::bind(&Phasefield<dim>::assemble_system_tangent_one_cell
                             ,this
                             ,param
                             ,std::placeholders::_1
                             ,std::placeholders::_2
                             ,std::placeholders::_3),
                   std::bind(&Phasefield<dim>::copy_local_to_global_K,
                             this,
                             std::placeholders::_1),
		   scratch_data,
                   per_task_data);

  timer.leave_subsection();
} 

template <int dim>
void Phasefield<dim>::copy_local_to_global_K(const PerTaskData_K &data)
{
  const unsigned int   dofs_per_cell = fe_m.dofs_per_cell;

  for (unsigned int i = 0; i < dofs_per_cell; ++i)
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
        tangent_matrix_m.add(data.local_dof_indices[i],
                             data.local_dof_indices[j],
                             data.cell_matrix(i, j));
}	

template <int dim>
void Phasefield<dim>::assemble_system_tangent_one_cell (const parameters::AllParameters &param,
                      		                        const typename DoFHandler<dim>::active_cell_iterator &cell
                                		       ,ScratchData_K &scratch
                                             	       ,PerTaskData_K &data)const
{
  data.reset();
  scratch.reset();
  scratch.fe_values.reinit(cell);
  
  const unsigned int   dofs_per_cell = fe_m.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula_m.size();

  double History;

  PointHistory *local_quadrature_point_data =
            reinterpret_cast<PointHistory *>(cell->user_pointer());

  scratch.fe_values[u_extractor].get_function_symmetric_gradients(scratch.update,scratch.epsilon_vals);
  scratch.fe_values[u_extractor].get_function_symmetric_gradients(scratch.old_solution_m,scratch.old_epsilon_vals);
  scratch.fe_values[d_extractor].get_function_values(scratch.update,scratch.d_vals);


  for (unsigned int q = 0; q < n_q_points; ++q)
  {

	if(current_time_m == param.time.delta_t)
		History = 0;
	else
		History = std::max(local_quadrature_point_data[q].history,get_history(param.materialmodel.lambda,param.materialmodel.mu,scratch.old_epsilon_vals[q]));

        local_quadrature_point_data[q].history = History;

	scratch.sigma_plus= get_stress_plus(param.materialmodel.lambda,param.materialmodel.mu,scratch.epsilon_vals[q]);
	scratch.BigC_plus = get_BigC_plus(param.materialmodel.lambda,param.materialmodel.mu,scratch.epsilon_vals[q]);
	scratch.BigC_minus = get_BigC_minus(param.materialmodel.lambda,param.materialmodel.mu,scratch.epsilon_vals[q]);

	for(unsigned int k=0;k<dofs_per_cell;++k)
	{
		scratch.shape_d[k] = scratch.fe_values[d_extractor].value(k,q);
		scratch.grad_shape_u[k] =  scratch.fe_values[u_extractor].gradient(k, q);
		scratch.grad_shape_d[k] =  scratch.fe_values[d_extractor].gradient(k, q);
		scratch.sym_grad_shape_u[k] =  scratch.fe_values[u_extractor].symmetric_gradient(k, q);
	}
            for (unsigned int i = 0; i < dofs_per_cell; ++i){

    		    for (unsigned int j = i; j < dofs_per_cell; ++j){

			if((dof_block_identifier_m[i] == d_dof_m) && (dof_block_identifier_m[j] == d_dof_m))
			{
                 	data.cell_matrix(i, j) +=( (param.pf.g_c/param.pf.l)*scratch.shape_d[i]*scratch.shape_d[j] 
						  + (param.pf.g_c*param.pf.l)* scalar_product(scratch.grad_shape_d[i],scratch.grad_shape_d[j])
					           + 2*scratch.shape_d[i]*scratch.shape_d[j]*History
					           + (param.materialmodel.viscosity/param.time.delta_t)*scratch.shape_d[i]*scratch.shape_d[j]) *scratch.fe_values.JxW(q);
			}
			else if((dof_block_identifier_m[i] == d_dof_m) && (dof_block_identifier_m[j] == u_dof_m))
			{
			data.cell_matrix(i,j) += 0;
			}
			else if((dof_block_identifier_m[i] == u_dof_m) && (dof_block_identifier_m[j] == d_dof_m))
			{
			data.cell_matrix(i,j) += ( 2*scalar_product(scratch.grad_shape_u[i],scratch.sigma_plus)*(1-scratch.d_vals[q])*scratch.shape_d[j])*scratch.fe_values.JxW(q);
			}
			else if((dof_block_identifier_m[i] == u_dof_m) && (dof_block_identifier_m[j] == u_dof_m))
			{
			double_contract(scratch.tmp1, scratch.BigC_plus, scratch.sym_grad_shape_u[j]);
			double_contract(scratch.tmp2, scratch.BigC_minus,scratch.sym_grad_shape_u[j]);

			data.cell_matrix(i,j) += (- scalar_product(scratch.grad_shape_u[i],scratch.tmp1)*(1 - scratch.d_vals[q])*(1 - scratch.d_vals[q])
					          - param.pf.k*scalar_product(scratch.grad_shape_u[i],scratch.tmp1)
					          - scalar_product(scratch.grad_shape_u[i],scratch.tmp2)  )*scratch.fe_values.JxW(q);
			}
                }
            }
          }

        for(unsigned int i=0;i<dofs_per_cell;++i)
            for(unsigned int j=0;j<i;++j)
                data.cell_matrix(i,j) = data.cell_matrix(j,i);

  cell->get_dof_indices(data.local_dof_indices);

}

template <int dim>
struct Phasefield<dim>::PerTaskData_rhs
{
    Vector<double>        cell_rhs;
    std::vector<types::global_dof_index> local_dof_indices;
    PerTaskData_rhs(const FiniteElement<dim> &fe)
      :
      cell_rhs(fe.dofs_per_cell),
      local_dof_indices(fe.dofs_per_cell)
    {}

    void reset()
    {
      cell_rhs = 0.0;
    }
};

template <int dim>
struct Phasefield<dim>::ScratchData_rhs
{
    FEValues<dim> fe_values;

    std::vector<SymmetricTensor<2,dim>> epsilon_vals;
    std::vector<SymmetricTensor<2,dim>> old_epsilon_vals;
    std::vector<double> d_vals;
    std::vector<double> old_d_vals;
    std::vector<Tensor<1,dim>> grad_d; 
    SymmetricTensor<2,dim> sigma_plus;
    SymmetricTensor<2,dim> sigma_minus;

    const BlockVector<double> & update;
    const BlockVector<double> & old_solution_m;
    std::vector<double> shape_d;
    std::vector<Tensor<1, dim>> shape_u;
    std::vector<Tensor<2, dim>> grad_shape_u;
    std::vector<Tensor<1, dim>> grad_shape_d;
   
    ScratchData_rhs( const FiniteElement<dim> &fe_cell
                  ,const QGauss<dim> &qf_cell
                  ,const UpdateFlags uf_cell
	          ,const BlockVector<double> &upd
		  ,const BlockVector<double> &old_sol_m)
     	        :fe_values(fe_cell, qf_cell, uf_cell)
		,epsilon_vals(qf_cell.size())
		,old_epsilon_vals(qf_cell.size())
		,d_vals(qf_cell.size())
		,old_d_vals(qf_cell.size())
	 	,grad_d(qf_cell.size())
		,update(upd)
		,old_solution_m(old_sol_m)
		,shape_d(fe_cell.dofs_per_cell)
		,shape_u(fe_cell.dofs_per_cell)
		,grad_shape_u(fe_cell.dofs_per_cell)
		,grad_shape_d(fe_cell.dofs_per_cell)
		{}

    ScratchData_rhs(const ScratchData_rhs &Scratch)
      		:fe_values(Scratch.fe_values.get_fe()
                ,Scratch.fe_values.get_quadrature()
                ,Scratch.fe_values.get_update_flags())
      		,epsilon_vals(Scratch.epsilon_vals)
		,old_epsilon_vals(Scratch.old_epsilon_vals)
      		,d_vals(Scratch.d_vals)
		,old_d_vals(Scratch.old_d_vals)
		,grad_d(Scratch.grad_d)
		,sigma_plus(Scratch.sigma_plus)
		,sigma_minus(Scratch.sigma_minus)
		,update(Scratch.update)
		,old_solution_m(Scratch.old_solution_m)
		,shape_d(Scratch.shape_d)
		,shape_u(Scratch.shape_u)
		,grad_shape_u(Scratch.grad_shape_u)
       		,grad_shape_d(Scratch.grad_shape_d)
		{}

    void reset()
    {
    const unsigned int n_q_points = epsilon_vals.size();
    for(unsigned int q=0;q<n_q_points;++q)
    {
	epsilon_vals[q] = 0.0;
    	old_epsilon_vals[q] = 0.0;
    	d_vals[q] = 0.0;
	old_d_vals[q]=0.0;
	grad_d[q]=0.0;
    }

    }
};

template <int dim>
void Phasefield<dim>::assemble_system_rhs (const AllParameters &param,BlockVector<double> & update)
{
  timer.enter_subsection("Assemble system rhs");

  system_rhs_m = 0.0;

  const UpdateFlags uf_cell(update_values    |  update_gradients |
                            update_quadrature_points |
			      update_JxW_values);

  PerTaskData_rhs per_task_data(fe_m);
  ScratchData_rhs scratch_data(fe_m, quadrature_formula_m, uf_cell,update,old_solution_m);

  WorkStream::run (dof_handler_m.begin_active(),
                   dof_handler_m.end(),
                   std::bind(&Phasefield<dim>::assemble_system_rhs_one_cell
                             ,this
                             ,param
                             ,std::placeholders::_1
                             ,std::placeholders::_2
                             ,std::placeholders::_3),
                   std::bind(&Phasefield<dim>::copy_local_to_global_rhs,
                             this,
                             std::placeholders::_1),
		   scratch_data,
                   per_task_data);

  timer.leave_subsection();
} 

template <int dim>
void Phasefield<dim>::copy_local_to_global_rhs(const PerTaskData_rhs &data)
{
  const unsigned int   dofs_per_cell = fe_m.dofs_per_cell;

  for (unsigned int i = 0; i < dofs_per_cell; ++i)
	  system_rhs_m(data.local_dof_indices[i]) +=  data.cell_rhs(i);
}	


template <int dim>
void Phasefield<dim>::assemble_system_rhs_one_cell (const parameters::AllParameters &param
                      		                    ,const typename DoFHandler<dim>::active_cell_iterator &cell
                                		    ,ScratchData_rhs &scratch
                                             	    ,PerTaskData_rhs &data)const
{
  data.reset();
  scratch.reset();
  scratch.fe_values.reinit(cell);
  cell->get_dof_indices(data.local_dof_indices);
  
  const unsigned int   dofs_per_cell = fe_m.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula_m.size();
  
  Tensor<1,dim> body_force;
  body_force = 0;
  
  double History;

  PointHistory *local_quadrature_point_data =
            reinterpret_cast<PointHistory *>(cell->user_pointer());

  scratch.fe_values[u_extractor].get_function_symmetric_gradients(scratch.update,scratch.epsilon_vals);
  scratch.fe_values[u_extractor].get_function_symmetric_gradients(scratch.old_solution_m,scratch.old_epsilon_vals);
  scratch.fe_values[d_extractor].get_function_gradients(scratch.update,scratch.grad_d);
  scratch.fe_values[d_extractor].get_function_values(scratch.update,scratch.d_vals);
  scratch.fe_values[d_extractor].get_function_values(scratch.old_solution_m,scratch.old_d_vals);


  for (unsigned int q = 0; q < n_q_points; ++q)
  {

	if(current_time_m == param.time.delta_t)
		History = 0;
	else
		History = std::max(local_quadrature_point_data[q].history,get_history(param.materialmodel.lambda,param.materialmodel.mu,scratch.old_epsilon_vals[q]));

        local_quadrature_point_data[q].history = History;

	scratch.sigma_plus= get_stress_plus(param.materialmodel.lambda,param.materialmodel.mu,scratch.epsilon_vals[q]);
	scratch.sigma_minus= get_stress_minus(param.materialmodel.lambda,param.materialmodel.mu,scratch.epsilon_vals[q]);


	for(unsigned int k=0;k<dofs_per_cell;++k)
	{
		scratch.shape_d[k] = scratch.fe_values[d_extractor].value(k,q);
		scratch.shape_u[k] = scratch.fe_values[u_extractor].value(k,q);
		scratch.grad_shape_u[k] =  scratch.fe_values[u_extractor].gradient(k, q);
		scratch.grad_shape_d[k] =  scratch.fe_values[d_extractor].gradient(k, q);

	}

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
	{
			if(dof_block_identifier_m[i] == d_dof_m )
                        {
                         data.cell_rhs(i) -= ( (param.pf.g_c/param.pf.l)*scratch.shape_d[i] *scratch.d_vals[q]
					+ (param.pf.g_c*param.pf.l)*scalar_product(scratch.grad_shape_d[i],scratch.grad_d[q])
					- 2*scratch.shape_d[i]*(1 - scratch.d_vals[q]) *History
					+ (param.materialmodel.viscosity/param.time.delta_t)*scratch.shape_d[i] * scratch.d_vals[q]
					- (param.materialmodel.viscosity/param.time.delta_t)*scratch.shape_d[i] * scratch.old_d_vals[q] 
					   )*scratch.fe_values.JxW(q);

                        }
                        else if(dof_block_identifier_m[i] == u_dof_m)
                        {
                        data.cell_rhs(i) -= (-scalar_product(scratch.grad_shape_u[i],scratch.sigma_plus) *(1 - scratch.d_vals[q])*(1-scratch.d_vals[q])
					- param.pf.k*scalar_product(scratch.grad_shape_u[i],scratch.sigma_plus)
					- scalar_product(scratch.grad_shape_u[i],scratch.sigma_minus)
					+ scalar_product(scratch.shape_u[i],body_force)
					 ) *scratch.fe_values.JxW(q);
                        }
    	}
  }

}

template <int dim>
void Phasefield<dim>::solve_nonlinear_newton(const AllParameters &param,
                                                 BlockVector<double> &solution_delta){

    BlockVector<double> newton_update(dofs_per_block_m);

    error_residual.reset();
    error_residual_0.reset();
    error_residual_norm.reset();
    
    error_update.reset();
    error_update_0.reset();
    error_update_norm.reset();


    print_header();
    unsigned int new_iter = 0;
    for (; new_iter < param.newtonraphson.max_new_ite; ++new_iter) {

        std::cout << " " << std::setw(2) << new_iter << " " << std::flush;
	
        tangent_matrix_m = 0.0;
        system_rhs_m = 0.0;

	if(current_time_m < param.time.time_change_interval)
	{
        make_constraints(new_iter,(param.time.delta_t/param.time.end_time),param.pf.u_total);

	}
	else
	{
	make_constraints(new_iter,(param.time.delta_t_f/param.time.end_time),param.pf.u_total);
	}

	BlockVector<double> current_sol = solution_m;//Includes previous timestep solution also
	current_sol += solution_delta;		//Gives updated solution till now 
	assemble_system_tangent(param,current_sol);	//Evaluation at solution calculated until now
      	assemble_system_rhs(param,current_sol);

	constraints_m.condense(tangent_matrix_m,system_rhs_m);		//as no distribute_l_to_g() in assemble_sys_tangent()

        get_error_residual(error_residual);

        if(new_iter==0)
            error_residual_0 = error_residual;

        error_residual_norm = error_residual;
        error_residual_norm.normalize(error_residual_0);
// checking
/*
	if(new_iter > 0 && error_residual_norm.u < param.newtonraphson.res_tol_u 
			&& error_residual_norm.d < param.newtonraphson.res_tol_d
			&& error_update_norm.u < param.newtonraphson.nu_tol_u
			&& error_update_norm.d < param.newtonraphson.nu_tol_d)
	{      	
            std::cout<<"Converged"<<std::endl;
            print_footer();
	    break;
        }
*/
//        assemble_system(param,current_sol/*solution_delta*/);		//Doesn't matter
//        make_constraints(new_iter,current_time_m);

        const std::pair<unsigned int,double>
                lin_solver_output = solve_linear_sys(param,newton_update);



	get_newton_update_error(newton_update,error_update);
	if(new_iter==0)
            error_update_0 = error_update;

        error_update_norm = error_update;
        error_update_norm.normalize(error_update_0);

        solution_delta += newton_update;

	
	//checking														
      	if(new_iter > 0 && (error_residual_norm.u < param.newtonraphson.res_tol_u || error_residual.u < param.newtonraphson.res_tol_u)
			&& (error_residual_norm.d < param.newtonraphson.res_tol_d || error_residual.d < param.newtonraphson.res_tol_d)
			&& (error_update_norm.u < param.newtonraphson.nu_tol_u || error_update.u < param.newtonraphson.nu_tol_u)
			&& (error_update_norm.d < param.newtonraphson.nu_tol_d || error_update.d < param.newtonraphson.nu_tol_d) )
	{      	
            std::cout<<"Converged"<<std::endl;
            print_footer();
	    break;
        }

        std::cout << " | " << std::fixed << std::setprecision(3) << std::setw(4)
                            << std::scientific << lin_solver_output.first << "   "
                            << lin_solver_output.second 
                            << "  " << error_residual.u  << "  " << error_residual_norm.u
			    << "  " << error_residual.d  << "  " << error_residual_norm.d
			    << "  " << error_update.u    << "  " << error_update_norm.u
			    << "  " << error_update.d    << "  " << error_update_norm.d
			    << "  " << get_energy(param,current_sol)
			    << "  " << std::endl;

    }
    AssertThrow (new_iter < param.newtonraphson.max_new_ite,
                   ExcMessage("No convergence in nonlinear solver!"));

}

template <int dim>
std::pair<unsigned int,double> Phasefield<dim>::solve_linear_sys(const AllParameters &param,BlockVector<double> &newton_update)
{
  timer.enter_subsection("Linear solver-d");
/*
  std::cout<<std::endl;
  std::cout<<"Initial Newton_update.block(d_dof_m): "<<newton_update.block(d_dof_m).l2_norm()<<std::endl;
  std::cout<<"Initial Newton_update.block(u_dof_m): "<<newton_update.block(u_dof_m).l2_norm()<<std::endl;
*/

  unsigned int lin_ite = 0;
  double lin_res = 0;

  Vector<double> tmp1(newton_update.block(u_dof_m).size());
  Vector<double> tmp2(newton_update.block(u_dof_m).size());

/*
  SolverControl           solver_control (tangent_matrix_m.block(d_dof_m,d_dof_m).m(),
                                          param.linearsolver.cg_tol*system_rhs_m.block(d_dof_m).l2_norm());
  GrowingVectorMemory<Vector<double> > GVM;
  SolverCG<Vector<double>>    cg (solver_control,GVM);
  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(tangent_matrix_m.block(d_dof_m,d_dof_m), param.linearsolver.relax_prm);
  cg.solve (tangent_matrix_m.block(d_dof_m,d_dof_m), newton_update.block(d_dof_m), system_rhs_m.block(d_dof_m),preconditioner);
*/

  SparseDirectUMFPACK k_dd;
  k_dd.initialize(tangent_matrix_m.block(d_dof_m,d_dof_m)); 
  k_dd.vmult(newton_update.block(d_dof_m),system_rhs_m.block(d_dof_m));
/*
  lin_ite = solver_control.last_step();
  lin_res = solver_control.last_value();
*/
 /*
  std::cout<<std::endl;
  std::cout<<"Solved Newton_update.block(d_dof_m)@1: "<<newton_update.block(d_dof_m).l2_norm()<<std::endl;
  std::cout<<"Initial Newton_update.block(u_dof_m)@1: "<<newton_update.block(u_dof_m).l2_norm()<<std::endl;
*/

 constraints_m.distribute (newton_update);		//?? //no effect
 timer.leave_subsection();
/*
  std::cout<<std::endl;
  std::cout<<"Solved Newton_update.block(d_dof_m)@2: "<<newton_update.block(d_dof_m).l2_norm()<<std::endl;
  std::cout<<"Initial Newton_update.block(u_dof_m)@2: "<<newton_update.block(u_dof_m).l2_norm()<<std::endl;
*/
  timer.enter_subsection("Linear solver-u");
  
  SparseDirectUMFPACK k_uu;
  k_uu.initialize(tangent_matrix_m.block(u_dof_m,u_dof_m)); 

//  std::cout<<"Error_d: tmp [f_d]: "<<system_rhs_m.block(d_dof_m).l2_norm()<<std::endl;

  tmp1 = system_rhs_m.block(u_dof_m);

//  std::cout<<"system_rhs_m.block(u_dof_m):tmp1 [f_u]: "<<tmp1.l2_norm()<<std::endl;

  tangent_matrix_m.block(u_dof_m,d_dof_m).vmult(tmp2,newton_update.block(d_dof_m));

//  std::cout<<"tmp2 :[k_ud*delta_d]: "<<tmp2.l2_norm()<<std::endl;

  tmp1 -= tmp2;

//  std::cout<<"tmp1: [f_u - k_ud*delta_d]: "<<tmp1.l2_norm()<<std::endl;

  k_uu.vmult(newton_update.block(u_dof_m),tmp1);

/*
  SolverControl           solver_control1 (tangent_matrix_m.block(u_dof_m,u_dof_m).m(),
                                          param.linearsolver.cg_tol*A.l2_norm());
  GrowingVectorMemory<Vector<double> > GVM1;
  SolverCG<Vector<double>>    cg1 (solver_control1,GVM1);
  PreconditionSSOR<> preconditioner1;
  preconditioner1.initialize(tangent_matrix_m.block(u_dof_m,u_dof_m), param.linearsolver.relax_prm);
  cg1.solve (tangent_matrix_m.block(u_dof_m,u_dof_m), newton_update.block(u_dof_m), A,preconditioner1);
*/
/*
  std::cout<<std::endl;
  std::cout<<"Solved Newton_update.block(d_dof_m): "<<newton_update.block(d_dof_m).l2_norm()<<std::endl;
  std::cout<<"Solved Newton_update.block(u_dof_m): "<<newton_update.block(u_dof_m).l2_norm()<<std::endl;
*/

  timer.leave_subsection();
  constraints_m.distribute (newton_update);
/*
  lin_ite = solver_control1.last_step();
  lin_res = solver_control1.last_value();
*/
/*
  std::cout<<std::endl;
  std::cout<<"Solved Newton_update.block(d_dof_m): "<<newton_update.block(d_dof_m).l2_norm()<<std::endl;
  std::cout<<"Solved Newton_update.block(u_dof_m): "<<newton_update.block(u_dof_m).l2_norm()<<std::endl;
*/
  return std::make_pair(lin_ite,lin_res);
}

template <int dim>
void Phasefield<dim>::output_results (const AllParameters &param,unsigned int cycle) const
{
    if((cycle%param.time.op_freq) == 0)
    {
    DataOut<dim> data_out;
    
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(dim,
                                  DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
    
    std::vector<std::string> solution_name(dim, "displacement");
    solution_name.emplace_back("phase_field");

    data_out.attach_dof_handler(dof_handler_m);
    data_out.add_data_vector(solution_m,
                             solution_name,
                             DataOut<dim>::type_dof_data,
                             data_component_interpretation);
    
    Vector<double> soln(solution_m.size());
    for (unsigned int i = 0; i < soln.size(); ++i)
      soln(i) = solution_m(i);
    
    MappingQEulerian<dim> q_mapping(param.fesys.fe_degree, dof_handler_m, soln);
    data_out.build_patches(q_mapping, param.fesys.fe_degree);
    std::ofstream output ("solution-"
                          + std::to_string(dim)
                          + "d-"
                          + std::to_string(cycle/param.time.op_freq)
                          + ".vtk");
    data_out.write_vtk(output);
    } 
}

template <int dim>
void Phasefield<dim>::import_mesh(const AllParameters &param){

    std::string grid_name;
    grid_name += param.geometrymodel.meshfile;

    GridIn<dim> grid_in;
    grid_in.attach_triangulation(triangulation_m);
    std::ifstream input_file(grid_name.c_str());
    Assert(dim==2, ExcInternalError());
//  grid_in.read_abaqus(input_file,false);
    grid_in.read_ucd(input_file,false);

    triangulation_m.refine_global (param.geometrymodel.gl_ref);
/*
    const bool write_grid = false;
    GridOut::OutputFormat meshOutputFormat = GridOut::vtk;
    if (write_grid)
    {
        const auto &      inputMeshFile =  param.geometrymodel.meshfile;
        GridOut  	  gridOut;
        const auto &      startPos = inputMeshFile.find_last_of("/\\") + 1;
        const auto &      endPos   = inputMeshFile.find_last_of('.');
        const std::string outMeshFileName =
                inputMeshFile.substr(startPos, endPos - startPos) +
                GridOut::default_suffix(meshOutputFormat);
        std::ofstream gridOutStream(outMeshFileName);
        gridOut.write(triangulation_m, gridOutStream, meshOutputFormat);
        std::cout << "The mesh has been written to " << outMeshFileName
                  << std::endl;
    }
*/

    set_boundary_id();
/*
    GridOut grid_out;

    GridOutFlags::Svg svg_flags;
    svg_flags.labeln_id = true;
    grid_out.set_flags(svg_flags);
    std::ofstream output_file(grid_name.c_str());
    grid_out.write_svg(triangulation_m,output_file);
*/
}     
 



template <int dim>
void Phasefield<dim>::print_header(){
    const unsigned int l_width = 130;
        for (unsigned int i = 0; i < l_width; ++i)
        {
            std::cout << "_";
        }
        std::cout << std::endl;

        std::cout << "SOLVER STEP"
                  << "| LIN_IT  LIN_RES     "
		  << "RES_U    RES_U_N    RES_D     RES_D_N     "                    
		  << "NU_U       NU_U_N     NU_D      NU_D_N      ENERGY"

		  << std::endl;

        for (unsigned int i = 0; i < l_width; ++i)
        {
            std::cout << "_";
        }
        std::cout << std::endl;

}

template <int dim>
void Phasefield<dim>::print_footer(){
    const unsigned int l_width = 130;
        for (unsigned int i = 0; i < l_width; ++i)
        {
            std::cout << "_";
        }
        std::cout << std::endl;

        std::cout << "Errors:" << std::endl
                  << "Res_u_n: \t\t" << error_residual_norm.u << std::endl
                  << "Res_d_n: \t\t" << error_residual_norm.d << std::endl
		  << std::endl;

}

template <int dim>
void Phasefield<dim>::get_error_residual(Error& error_residual){


    BlockVector<double> err_res(dofs_per_block_m);
    
    for (unsigned int i=0;i<dof_handler_m.n_dofs();++i) {
        if(!constraints_m.is_constrained(i)){
            err_res(i)=system_rhs_m(i);
           }

        }
    error_residual.norm = err_res.l2_norm();
    error_residual.u = err_res.block(u_dof_m).l2_norm();
    error_residual.d = err_res.block(d_dof_m).l2_norm();

}

template <int dim>
void Phasefield<dim>::get_newton_update_error(const BlockVector<double> &newton_update
					     ,Error& error_update){

    BlockVector<double> err_update(dofs_per_block_m);
    for (unsigned int i=0;i<dof_handler_m.n_dofs();++i) {
        if(!constraints_m.is_constrained(i)){
            err_update(i)= newton_update(i);
           }

        }
    error_update.norm = err_update.l2_norm();
    error_update.u = err_update.block(u_dof_m).l2_norm();
    error_update.d = err_update.block(d_dof_m).l2_norm();

}

template <int dim>
void Phasefield<dim>::project_back_phase_field ()
{
  typename DoFHandler<dim>::active_cell_iterator cell =
    dof_handler_m.begin_active(), endc = dof_handler_m.end();

//  int count1=0,count2=0;
  
  std::vector<unsigned int> local_dof_indices(fe_m.dofs_per_cell);
  for (; cell != endc; ++cell)
      {
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i=0; i<fe_m.dofs_per_cell; ++i)
          {
            const unsigned int comp_i = fe_m.system_to_component_index(i).first;
              
//	    const unsigned int k_group = fe_m.system_to_base_index(i).first.first;//Return for shape function index(i) the base element it belongs to.
	      
/*	    if (k_group == u_dof_m)
	    {
		   ++ count1;
		   continue;
	    }
*/
  	    if (comp_i != dim)
              continue; // only look at phase field

            const unsigned int idx = local_dof_indices[i];
	
//	    if(solution_m(idx) > 1.0 )
//	    std::cout<<"Before solution_m(idx): "<<solution_m(idx)<<std::endl;
            
	    solution_m(idx) = std::max(0.0,std::min(static_cast<double>(solution_m(idx)), 1.0));
//            ++ count2 ;

//	    if(solution_m(idx) > 1.0 )
//	    std::cout<<"After solution_m(idx): "<<solution_m(idx)<<std::endl;

	  }
      }
//  std::cout<<"count u:"<<count1<<std::endl;
//  std::cout<<"count d:"<<count2<<std::endl;
}


template <int dim>
void Phasefield<dim>::make_constraints(unsigned int &itr,const double load_ratio,const double u_total){
    std::cout<<" CST "<<std::flush;

    if(itr>1){				//no effect of removing
        return;
    }

   constraints_m.clear();
/*
   DoFTools::make_hanging_node_constraints (dof_handler_m,
                                           constraints_m);
*/
   std::vector<bool> component_mask(dim+1, false);
     
   //Tension test
   { 
   component_mask[0] = false;
   component_mask[1] = true;
   VectorTools::interpolate_boundary_values(dof_handler_m, 2,
                                               ZeroFunction<dim>(n_components_m), constraints_m, component_mask);
   }
   
   {
   component_mask[0] = true;
   component_mask[1] = true;
   VectorTools::interpolate_boundary_values(dof_handler_m, 3,
                                               BoundaryTension<dim>(itr,load_ratio,u_total), constraints_m, component_mask);
   }
 //?? u_x = traction free (false) 
 
  /*
   // Shear Test{
   { 
   component_mask[0] = false;
   component_mask[1] = true;
   VectorTools::interpolate_boundary_values(dof_handler_m, 0,
                                               ZeroFunction<dim>(n_components_m), constraints_m, component_mask);
   }
   { 
   component_mask[0] = false;
   component_mask[1] = true;
   VectorTools::interpolate_boundary_values(dof_handler_m, 1,
                                               ZeroFunction<dim>(n_components_m), constraints_m, component_mask);
   }
   { 
   component_mask[0] = true;
   component_mask[1] = true;
   VectorTools::interpolate_boundary_values(dof_handler_m, 2,
                                               ZeroFunction<dim>(n_components_m), constraints_m, component_mask);
   }
   {
   component_mask[0] = true;
   component_mask[1] = true;
   VectorTools::interpolate_boundary_values(dof_handler_m, 3,
                                               BoundaryShear<dim>(itr,load_ratio,u_total), constraints_m, component_mask);
   }
   // Lower part of slit
   {
   component_mask[0] = false;
   component_mask[1] = true;
   VectorTools::interpolate_boundary_values(dof_handler_m, 4,
                                               ZeroFunction<dim>(n_components_m), constraints_m, component_mask);
   }
   //shear end}
   */

   constraints_m.close ();

     

}

template <int dim>
void Phasefield<dim>::compute_lefm_errors(const AllParameters &param)
{
    const ComponentSelectFunction<dim> displacement_mask(std::make_pair(0, dim),dim + 1);
    
    Vector<float> difference_per_cell(triangulation_m.n_active_cells());
    
    VectorTools::integrate_difference(dof_handler_m,
                                      solution_m,
                                      Reference_solution<dim>(),
                                      difference_per_cell,
                                      QGauss<dim>(param.fesys.quad_order),
                                      VectorTools::L2_norm,
                                      & displacement_mask);
    
    const double Displacement_L2_error =
      VectorTools::compute_global_error(triangulation_m,
                                        difference_per_cell,
                                        VectorTools::L2_norm);

    std::cout << std::endl<< "   Displacement L2 Error: " << Displacement_L2_error << std::endl;
}


template <int dim>
void Phasefield<dim>::run(const AllParameters &param){


    using namespace constants;

    BlockVector<double>       solution_delta(dofs_per_block_m);

    long double present_time_tol;
    present_time_tol = time_tol_fac * param.time.delta_t;
    current_time_m = param.time.start_time;
    unsigned int current_timestep = 0;

    std::cout<<std::endl;
    std::cout<<"Min cell dia(before local ref):"<<GridTools::minimal_cell_diameter(triangulation_m)<<std::endl;

	
    //Local Pre-refinement for Tension test
    for(unsigned int i=0;i<param.geometrymodel.lc_ref;++i)
    {
    typename DoFHandler<dim>::active_cell_iterator cell =
        dof_handler_m.begin_active(), endc = dof_handler_m.end();

      for (; cell != endc; ++cell)
         {
            for (unsigned int vertex = 0;vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex)
              {
                Tensor<1, dim> cell_vertex = (cell->vertex(vertex));
                if (cell_vertex[0] <= 0.6 && cell_vertex[0] >= 0.0
                    && cell_vertex[1] <= 0.65 && cell_vertex[1] >= 0.35)
 		{
                    cell->set_refine_flag();
                    break;
                }
	      }
          }

     for (; cell != endc; ++cell)
        if (cell->level() == static_cast<int>(param.geometrymodel.gl_ref+4))
          cell->clear_refine_flag();


//    triangulation_m.prepare_coarsening_and_refinement();
    triangulation_m.execute_coarsening_and_refinement();//execute_refinement();
    setup_system();
    setup_quadrature_point_history();
    solution_delta.reinit(dofs_per_block_m);
    }
    
//    std::cout<<"Before initial crack:solution(d_dof_m).norm:"<<solution_m.block(d_dof_m).l2_norm()<<std::endl;
//    std::cout<<"Before initial crack:solution(u_dof_m).norm:"<<solution_m.block(u_dof_m).l2_norm()<<std::endl;

/*    //No effect/////////
    constraints_m.close();								 			//althoughalreadycalledinsetup
    VectorTools::interpolate(dof_handler_m
                         ,InitialCrack<dim>(GridTools::minimal_cell_diameter(triangulation_m)),old_solution_m); //??
    constraints_m.distribute(old_solution_m);						 //tomakeresultcontinuousagain
    solution_m = old_solution_m;
    ///////////////////

    std::cout<<"After initial crack: solution(d_dof_m).norm:"<<solution_m.block(d_dof_m).l2_norm()<<std::endl;
    std::cout<<"After initial crack: solution(u_dof_m).norm:"<<solution_m.block(u_dof_m).l2_norm()<<std::endl;
    std::cout<<std::endl;
    std::cout<<"After: solution(d_dof_m).norm.linfinity:"<<solution_m.block(d_dof_m).linfty_norm()<<std::endl;
*/

    output_results(param,current_timestep);
    current_time_m += param.time.delta_t;
    current_timestep++;

/////////Printing////////
	std::cout<<std::endl;
	std::cout<<"Min cell dia(after local ref):"<<GridTools::minimal_cell_diameter(triangulation_m)<<std::endl;
	std::cout<<"Parameter Global Ref: "<<param.geometrymodel.gl_ref<<std::endl;
	std::cout<<"Parameter Locl Ref: "<<param.geometrymodel.lc_ref<<std::endl;
	std::cout<<"Parameter vis: "<<param.materialmodel.viscosity<<std::endl;
	std::cout<<"Parameter g_c: "<<param.pf.g_c<<std::endl;
	std::cout<<"Parameter l: "<<param.pf.l<<std::endl;
	std::cout<<"Parameter k: "<<param.pf.k<<std::endl;
	std::cout<<"Parameter delta_t: "<<param.time.delta_t<<std::endl;
	std::cout<<"Parameter end_time: "<<param.time.end_time<<std::endl;
	std::cout<<"Parameter u_total: "<<param.pf.u_total<<std::endl;
	std::cout<<"Parameter time change interval: "<<param.time.time_change_interval<<std::endl;
	std::cout<<"Parameter Newton: res tol_u: "<<param.newtonraphson.res_tol_u<<std::endl;
	std::cout<<"Parameter Newton: res tol_d: "<<param.newtonraphson.res_tol_d<<std::endl;
	std::cout<<"Parameter Newton: new_upd tol_u: "<<param.newtonraphson.nu_tol_u<<std::endl;
	std::cout<<"Parameter Newton: new_upd tol_d: "<<param.newtonraphson.nu_tol_d<<std::endl;
	std::cout<<std::endl;
///////////////////////
   int counter=0; 
    
    while (current_time_m < param.time.end_time + present_time_tol)
	  
    {
	std::cout<<"Current time step: "<<current_timestep<<std::endl;

//	constraints_m.distribute(old_solution_m);

	solution_delta = 0.0;
        solve_nonlinear_newton(param,solution_delta);
        solution_m += solution_delta;	

//	project_back_phase_field();

//	compute_lefm_errors(param);			//For Reference solution

	output_results(param,current_timestep);

       	old_solution_m = solution_m;

	if((current_timestep%param.time.op_freq) == 0)
	{
	    statistics.add_value("Time",current_time_m);
	    compute_load(param.materialmodel.lambda,param.materialmodel.mu,solution_m);
       
	    std::ofstream stat_file ("statistics_m");
            statistics.write_text (stat_file,
                                   TableHandler::simple_table_with_separate_column_description);
            stat_file.close();
	counter++;
//	std::cout<<"Current time:"<<current_time_m<<"  counter:"<<counter<<std::endl;
//	std::cout<<"solution(d_dof_m).norm:"<<solution_m.block(d_dof_m).l2_norm()<<std::endl;
//	std::cout<<"solution(u_dof_m).norm:"<<solution_m.block(u_dof_m).l2_norm()<<std::endl;
//	std::cout<<std::endl;
	}
	if(current_time_m < param.time.time_change_interval)
	{
		current_time_m += param.time.delta_t;
	}
	else
	{
		current_time_m += param.time.delta_t_f;
	}
  	++current_timestep;
   }

}

template class thesis::Phasefield<2>;
template class thesis::Phasefield<3>;

