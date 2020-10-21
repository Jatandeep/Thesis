
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
  , timer(std::cout,TimerOutput::summary,TimerOutput::wall_times)
				
{
    import_mesh(param);
    dof_handler_m.initialize(triangulation_m, fe_m);
    setup_system();			
    setup_qph();
    determine_comp_extractor();
}



template <int dim>
Phasefield<dim>::~Phasefield ()
{
  dof_handler_m.clear ();
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

  constraints_m.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler_m,
                                           constraints_m);
  //get_constrained_initial_d();
    
  constraints_m.close ();

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
void Phasefield<dim>::set_boundary_id(const AllParameters &param)
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
              if(param.mod_strategy.comp_strategy=="normal")
              {
		            //left boundary
                if ((face_center[0] < 0.0+tol_machine) && (face_center[0] > 0.0-tol_machine)
                   )
                  cell->face(f)->set_boundary_id(1);
		            //right boundary
                else if ((face_center[0] < 1.0+tol_machine) && (face_center[0] > 1.0-tol_machine)
                        )
                  cell->face(f)->set_boundary_id(2);
		            // bottom boundary
                else if ((face_center[1] < 0.0+tol_machine) && (face_center[1] > 0.0-tol_machine)
                        )
                  cell->face(f)->set_boundary_id(3);
              	// top boundary
                else if ((face_center[1] < 1.0+tol_machine) && (face_center[1] > 1.0-tol_machine)
                        )
                  cell->face(f)->set_boundary_id(4);
              }
              if(param.mod_strategy.comp_strategy=="lefm")
              {
                //left boundary
                if ((face_center[0] < -0.5+tol_machine) && (face_center[0] > -0.5-tol_machine)
                   )
                  cell->face(f)->set_boundary_id(1);
		            //right boundary
                else if ((face_center[0] < 0.5+tol_machine) && (face_center[0] > 0.5-tol_machine)
                        )
                  cell->face(f)->set_boundary_id(2);
		            // bottom boundary
                else if ((face_center[1] < -0.5+tol_machine) && (face_center[1] > -0.5-tol_machine)
                        )
                  cell->face(f)->set_boundary_id(3);
              	// top boundary
                else if ((face_center[1] < 0.5+tol_machine) && (face_center[1] > 0.5-tol_machine)
                        )
                  cell->face(f)->set_boundary_id(4);
              }   
            }
          }

}

template<int dim>
void Phasefield<dim>::setup_qph()
{
  std::cout<<"CellDataStorage:setting up quadrature point data.."<<std::endl;

  const unsigned int   n_q_points    = quadrature_formula_m.size();
  quadrature_point_history.initialize(triangulation_m.begin_active()
                                      ,triangulation_m.end()
                                      ,n_q_points);
 
  typename Triangulation<dim>::active_cell_iterator cell =
      triangulation_m.begin_active(), endc = triangulation_m.end();

  for (;cell != endc; ++cell)
  {
      std::vector<std::shared_ptr<PointHistory>> lqph = quadrature_point_history.get_data(cell);
      Assert(lqph.size() == n_q_points, ExcInternalError());
  }

}

template <int dim>
void Phasefield<dim>::get_energy_v(const AllParameters &param
                                                     ,BlockVector<double> & update)
{
  FEValues<dim> fe_values (fe_m, quadrature_formula_m,
                           update_values   | update_gradients |
                           update_quadrature_points | update_JxW_values);

  const unsigned int   n_q_points    = quadrature_formula_m.size();

  double Elastic_energy=0,Fracture_energy=0,Total_energy=0;
    
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
	
		  Elastic_energy += ( (get_deg_func(d_vals[q]) + param.pf.k)*get_energy_density_plus(param.materialmodel.lambda,param.materialmodel.mu,epsilon_vals[q]) 
				                + get_energy_density_minus(param.materialmodel.lambda,param.materialmodel.mu,epsilon_vals[q]) )
                        *fe_values.JxW(q);

      Fracture_energy += ( param.pf.g_c*( (0.5/param.pf.l)*d_vals[q]*d_vals[q] 
                          + 0.5*param.pf.l* (scalar_product(grad_d[q],grad_d[q])) ))
                          *fe_values.JxW(q);
      }
    }
  Total_energy = Elastic_energy + Fracture_energy;

  statistics.add_value("Elastic Energy",Elastic_energy);
  statistics.add_value("Fracture Energy",Fracture_energy);
  statistics.add_value("Total Energy",Total_energy);
}

template <int dim>
std::pair<double,double> Phasefield<dim>::get_energy_p(const AllParameters &param
                                                     ,BlockVector<double> & update)
{
  FEValues<dim> fe_values (fe_m, quadrature_formula_m,
                           update_values   | update_gradients |
                           update_quadrature_points | update_JxW_values);

  const unsigned int   n_q_points    = quadrature_formula_m.size();

  double Elastic_energy=0,Fracture_energy=0,Total_energy=0;
    
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
	
		  Elastic_energy += ( (get_deg_func(d_vals[q]) + param.pf.k)*get_energy_density_plus(param.materialmodel.lambda,param.materialmodel.mu,epsilon_vals[q]) 
				                + get_energy_density_minus(param.materialmodel.lambda,param.materialmodel.mu,epsilon_vals[q]) )
                        *fe_values.JxW(q);

      Fracture_energy += ( param.pf.g_c*( (0.5/param.pf.l)*d_vals[q]*d_vals[q] 
                          + 0.5*param.pf.l* (scalar_product(grad_d[q],grad_d[q])) ))
                          *fe_values.JxW(q);
      }
    }
  Total_energy = Elastic_energy + Fracture_energy;

return std::make_pair(Total_energy,Fracture_energy);
}


template <int dim>
struct Phasefield<dim>::PerTaskData_d
{
    FullMatrix<double>        cell_matrix;
    Vector<double>        cell_rhs;

    std::vector<types::global_dof_index> local_dof_indices;
    PerTaskData_d(const FiniteElement<dim> &fe)
      :
      cell_matrix(fe.dofs_per_cell, fe.dofs_per_cell)
      ,cell_rhs(fe.dofs_per_cell)
      ,local_dof_indices(fe.dofs_per_cell)
    {}

    void reset()
    {
      cell_matrix = 0.0;
      cell_rhs = 0.0;
    }
};

template <int dim>
struct Phasefield<dim>::ScratchData_d
{
    FEValues<dim> fe_values;

    std::vector<SymmetricTensor<2,dim>> old_epsilon_vals;
    std::vector<double> d_vals;
    std::vector<Tensor<1,dim>> grad_d;
    std::vector<double> old_d_vals;
    
    const BlockVector<double> & update;
    const BlockVector<double> & old_solution_m;
    const double delta_t;
    std::vector<double> shape_d;
    std::vector<Tensor<1, dim>> grad_shape_d;

    ScratchData_d( const FiniteElement<dim> &fe_cell
                  ,const QGauss<dim> &qf_cell
                  ,const UpdateFlags uf_cell
	                ,const BlockVector<double> &upd
		              ,const BlockVector<double> &old_sol_m
                  ,const double delta_t)
     	            :fe_values(fe_cell, qf_cell, uf_cell)
                  ,old_epsilon_vals(qf_cell.size())
                  ,d_vals(qf_cell.size())
                  ,grad_d(qf_cell.size())
                  ,old_d_vals(qf_cell.size())
                  ,update(upd)
                  ,old_solution_m(old_sol_m)
                  ,delta_t(delta_t)
                  ,shape_d(fe_cell.dofs_per_cell)
                  ,grad_shape_d(fe_cell.dofs_per_cell)
                  {}

    ScratchData_d(const ScratchData_d &Scratch)
      		      :fe_values(Scratch.fe_values.get_fe()
                ,Scratch.fe_values.get_quadrature()
                ,Scratch.fe_values.get_update_flags())
      		  		,old_epsilon_vals(Scratch.old_epsilon_vals)
      		      ,d_vals(Scratch.d_vals)
                ,grad_d(Scratch.grad_d)
                ,old_d_vals(Scratch.old_d_vals)
                ,update(Scratch.update)
                ,old_solution_m(Scratch.old_solution_m)
                ,delta_t(Scratch.delta_t)
                ,shape_d(Scratch.shape_d)
                ,grad_shape_d(Scratch.grad_shape_d)
            		{}

    void reset()
    {
    const unsigned int n_q_points = d_vals.size();
    for(unsigned int q=0;q<n_q_points;++q)
    {
    	old_epsilon_vals[q] = 0.0;
    	d_vals[q] = 0.0;
      grad_d[q] =0.0;
      old_d_vals[q] = 0.0;
    }

    }
};

template <int dim>
void Phasefield<dim>::assemble_system_d (const AllParameters &param
                                        ,BlockVector<double> & update
                                        ,double delta_t)
{
  timer.enter_subsection("Assemble d");

  tangent_matrix_m = 0.0;
  system_rhs_m = 0.0;

  const UpdateFlags uf_cell(update_values    |  update_gradients |
                            update_quadrature_points |
			                      update_JxW_values);

  PerTaskData_d per_task_data(fe_m);
  ScratchData_d scratch_data(fe_m, quadrature_formula_m, uf_cell,update,old_solution_m,delta_t);

  WorkStream::run (dof_handler_m.begin_active(),
                   dof_handler_m.end(),
                   std::bind(&Phasefield<dim>::assemble_system_d_one_cell
                             ,this
                             ,param
                             ,std::placeholders::_1
                             ,std::placeholders::_2
                             ,std::placeholders::_3),
                   std::bind(&Phasefield<dim>::copy_local_to_global_d,
                             this,
                             std::placeholders::_1),
		               scratch_data,
                   per_task_data);

  timer.leave_subsection();
} 

template <int dim>
void Phasefield<dim>::copy_local_to_global_d(const PerTaskData_d &data)
{
  constraints_m.distribute_local_to_global(data.cell_matrix,data.cell_rhs,data.local_dof_indices
                                          ,tangent_matrix_m,system_rhs_m);	
}

template <int dim>
void Phasefield<dim>::assemble_system_d_one_cell (const parameters::AllParameters &param
                      		                        ,const typename DoFHandler<dim>::active_cell_iterator &cell
                                		              ,ScratchData_d &scratch
                                             	    ,PerTaskData_d &data)const
{
  data.reset();
  scratch.reset();
  scratch.fe_values.reinit(cell);
  
  const unsigned int   dofs_per_cell = fe_m.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula_m.size();

  double History;

  std::vector<std::shared_ptr<const PointHistory>> lqph = quadrature_point_history.get_data(cell);
  
  scratch.fe_values[u_extractor].get_function_symmetric_gradients(scratch.old_solution_m,scratch.old_epsilon_vals);
  scratch.fe_values[d_extractor].get_function_values(scratch.update,scratch.d_vals);
  scratch.fe_values[d_extractor].get_function_gradients(scratch.update,scratch.grad_d);
  scratch.fe_values[d_extractor].get_function_values(scratch.old_solution_m,scratch.old_d_vals);


  for (unsigned int q = 0; q < n_q_points; ++q)
  {

	if(current_time_m == param.time.delta_t)
		History = 0;
	else
		History = std::max(lqph[q]->history
                      ,get_history(param.materialmodel.lambda,param.materialmodel.mu,scratch.old_epsilon_vals[q]));

      //  if(cell==dof_handler_m.begin_active() && q==2){
      //   std::cout<<std::endl;
      //   std::cout<<"vertex_index: "<<cell->vertex_index(q)<<std::endl;
      //   std::cout<<"get_history(): "<<get_history(param.materialmodel.lambda,param.materialmodel.mu,scratch.old_epsilon_vals[q])<<std::endl;
      //   //std::cout<<"local_quadrature_point_data[q].history "<<local_quadrature_point_data[q].history<<std::endl;
      //   std::cout<<"lqph[q]->history:"<<lqph[q]->history<<std::endl;
      //   std::cout<<"History: "<<History<<std::endl;        
      //      std::cout<<"delta t inside d assembly: "<<scratch.delta_t<<std::endl;
      //std::cout<<"viscosity: "<<param.materialmodel.viscosity<<std::endl;
      //         }      

    lqph[q]->history = History; 

    for(unsigned int k=0;k<dofs_per_cell;++k)
    {
      scratch.shape_d[k] = scratch.fe_values[d_extractor].value(k,q);
      scratch.grad_shape_d[k] =  scratch.fe_values[d_extractor].gradient(k, q);
    }
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
    		    for (unsigned int j = i; j < dofs_per_cell; ++j)
            {

              if((dof_block_identifier_m[i] == d_dof_m) && (dof_block_identifier_m[j] == d_dof_m))
              {
                  data.cell_matrix(i, j) +=( (param.pf.g_c/param.pf.l)*scratch.shape_d[i]*scratch.shape_d[j] 
                                          + (param.pf.g_c*param.pf.l)* scalar_product(scratch.grad_shape_d[i],scratch.grad_shape_d[j])
                                          + 2*scratch.shape_d[i]*scratch.shape_d[j]*History
                                          + (param.materialmodel.viscosity/scratch.delta_t)*scratch.shape_d[i]*scratch.shape_d[j])
                                          *scratch.fe_values.JxW(q);
              }
              else if((dof_block_identifier_m[i] == d_dof_m) && (dof_block_identifier_m[j] == u_dof_m))
              {
                  data.cell_matrix(i,j) += 0;
              }
		
            }
            
            if(dof_block_identifier_m[i] == d_dof_m )
            {
                 data.cell_rhs(i) -= ( (param.pf.g_c/param.pf.l)*scratch.shape_d[i] *scratch.d_vals[q]
                                      + (param.pf.g_c*param.pf.l)*scalar_product(scratch.grad_shape_d[i],scratch.grad_d[q])
                                      - 2*scratch.shape_d[i]*(1 - scratch.d_vals[q]) *History
                                      + (param.materialmodel.viscosity/scratch.delta_t)*scratch.shape_d[i] * scratch.d_vals[q]
                                      - (param.materialmodel.viscosity/scratch.delta_t)*scratch.shape_d[i] * scratch.old_d_vals[q] 
                                        )*scratch.fe_values.JxW(q);

            }
      }
  }

  for(unsigned int i=0;i<dofs_per_cell;++i)
      for(unsigned int j=0;j<i;++j)
          data.cell_matrix(i,j) = data.cell_matrix(j,i);

  cell->get_dof_indices(data.local_dof_indices);

}

template <int dim>
struct Phasefield<dim>::PerTaskData_u
{
    FullMatrix<double>        cell_matrix;
    Vector<double>        cell_rhs;
    std::vector<types::global_dof_index> local_dof_indices;
    PerTaskData_u(const FiniteElement<dim> &fe)
      :
      cell_matrix(fe.dofs_per_cell,fe.dofs_per_cell)
      ,cell_rhs(fe.dofs_per_cell)
      ,local_dof_indices(fe.dofs_per_cell)
    {}

    void reset()
    {
      cell_matrix = 0.0;
      cell_rhs = 0.0;
    }
};

template <int dim>
struct Phasefield<dim>::ScratchData_u
{
    FEValues<dim> fe_values;

    std::vector<SymmetricTensor<2,dim>> epsilon_vals;
    std::vector<double> d_vals;
    SymmetricTensor<2,dim> sigma_plus;
    SymmetricTensor<2,dim> sigma_minus;
    SymmetricTensor<2,dim> tmp1,tmp2;
    SymmetricTensor<4,dim> BigC_plus;
    SymmetricTensor<4,dim> BigC_minus;

    const BlockVector<double> & update;
    std::vector<Tensor<1, dim>> shape_u;
    std::vector<Tensor<2, dim>> grad_shape_u;
    std::vector<SymmetricTensor<2, dim>> sym_grad_shape_u;
   
    ScratchData_u( const FiniteElement<dim> &fe_cell
                  ,const QGauss<dim> &qf_cell
                  ,const UpdateFlags uf_cell
	                ,const BlockVector<double> &upd)
     	            :fe_values(fe_cell, qf_cell, uf_cell)
                  ,epsilon_vals(qf_cell.size())
                  ,d_vals(qf_cell.size())
                  ,update(upd)
                  ,shape_u(fe_cell.dofs_per_cell)
                  ,grad_shape_u(fe_cell.dofs_per_cell)
                  ,sym_grad_shape_u(fe_cell.dofs_per_cell)
                  {}

    ScratchData_u(const ScratchData_u &Scratch)
      		        :fe_values(Scratch.fe_values.get_fe()
                  ,Scratch.fe_values.get_quadrature()
                  ,Scratch.fe_values.get_update_flags())
                  ,epsilon_vals(Scratch.epsilon_vals)
                  ,d_vals(Scratch.d_vals)
                  ,sigma_plus(Scratch.sigma_plus)
                  ,sigma_minus(Scratch.sigma_minus)
                  ,tmp1(Scratch.tmp1)
                  ,tmp2(Scratch.tmp2)
                  ,BigC_plus(Scratch.BigC_plus)
		              ,BigC_minus(Scratch.BigC_minus)
                  ,update(Scratch.update)
                  ,shape_u(Scratch.shape_u)
                  ,grad_shape_u(Scratch.grad_shape_u)
                  ,sym_grad_shape_u(Scratch.sym_grad_shape_u)
                  {}

    void reset()
    {
    const unsigned int n_q_points = epsilon_vals.size();
    for(unsigned int q=0;q<n_q_points;++q)
    {
	    epsilon_vals[q] = 0.0;
    	d_vals[q] = 0.0;
	  }

    }
};

template <int dim>
void Phasefield<dim>::assemble_system_u (const AllParameters &param,BlockVector<double> & update)
{
  timer.enter_subsection("Assemble u");

  tangent_matrix_m = 0.0;
  system_rhs_m = 0.0;

  const UpdateFlags uf_cell(update_values    |  update_gradients |
                            update_quadrature_points |
			                      update_JxW_values);

  PerTaskData_u per_task_data(fe_m);
  ScratchData_u scratch_data(fe_m, quadrature_formula_m, uf_cell,update);

  WorkStream::run (dof_handler_m.begin_active(),
                   dof_handler_m.end(),
                   std::bind(&Phasefield<dim>::assemble_system_u_one_cell
                             ,this
                             ,param
                             ,std::placeholders::_1
                             ,std::placeholders::_2
                             ,std::placeholders::_3),
                   std::bind(&Phasefield<dim>::copy_local_to_global_u,
                             this,
                             std::placeholders::_1),
		               scratch_data,
                   per_task_data);

  timer.leave_subsection();
} 

template <int dim>
void Phasefield<dim>::copy_local_to_global_u(const PerTaskData_u &data)
{
//  const unsigned int   dofs_per_cell = fe_m.dofs_per_cell;

 constraints_m.distribute_local_to_global(data.cell_matrix,data.cell_rhs,data.local_dof_indices
                                          ,tangent_matrix_m,system_rhs_m);	
}

template <int dim>
void Phasefield<dim>::assemble_system_u_one_cell (const parameters::AllParameters &param
                      		                        ,const typename DoFHandler<dim>::active_cell_iterator &cell
                                		              ,ScratchData_u &scratch
                                             	    ,PerTaskData_u &data)const
{
  data.reset();
  scratch.reset();
  scratch.fe_values.reinit(cell);
    
  const unsigned int   dofs_per_cell = fe_m.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula_m.size();
  
  Tensor<1,dim> body_force;
  body_force = 0;
  
  scratch.fe_values[u_extractor].get_function_symmetric_gradients(scratch.update,scratch.epsilon_vals);
  scratch.fe_values[d_extractor].get_function_values(scratch.update,scratch.d_vals);

  for (unsigned int q = 0; q < n_q_points; ++q)
  {

	  scratch.sigma_plus= get_stress_plus(param.materialmodel.lambda,param.materialmodel.mu,scratch.epsilon_vals[q]);
    scratch.sigma_minus= get_stress_minus(param.materialmodel.lambda,param.materialmodel.mu,scratch.epsilon_vals[q]);
    scratch.BigC_plus = get_BigC_plus(param.materialmodel.lambda,param.materialmodel.mu,scratch.epsilon_vals[q]);
    scratch.BigC_minus = get_BigC_minus(param.materialmodel.lambda,param.materialmodel.mu,scratch.epsilon_vals[q]);

    for(unsigned int k=0;k<dofs_per_cell;++k)
    {
      scratch.shape_u[k] = scratch.fe_values[u_extractor].value(k,q);
      scratch.grad_shape_u[k] =  scratch.fe_values[u_extractor].gradient(k, q);
      scratch.sym_grad_shape_u[k] =  scratch.fe_values[u_extractor].symmetric_gradient(k, q);
    }

    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {

      for (unsigned int j = i; j < dofs_per_cell; ++j)
      {

        if((dof_block_identifier_m[i] == u_dof_m) && (dof_block_identifier_m[j] == d_dof_m))
        {
          data.cell_matrix(i,j) += 0.0;
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

			
      if(dof_block_identifier_m[i] == u_dof_m)
      {
          data.cell_rhs(i) -= (-scalar_product(scratch.grad_shape_u[i],scratch.sigma_plus) *(1 - scratch.d_vals[q])*(1-scratch.d_vals[q])
                              - param.pf.k*scalar_product(scratch.grad_shape_u[i],scratch.sigma_plus)
                              - scalar_product(scratch.grad_shape_u[i],scratch.sigma_minus)
                              + scalar_product(scratch.shape_u[i],body_force)
                              ) *scratch.fe_values.JxW(q);
      }
    }
  }
  for(unsigned int i=0;i<dofs_per_cell;++i)
      for(unsigned int j=0;j<i;++j)
          data.cell_matrix(i,j) = data.cell_matrix(j,i);

  cell->get_dof_indices(data.local_dof_indices);
}

template <int dim>
unsigned int Phasefield<dim>::solve_nonlinear_newton(const AllParameters &param
                                            ,BlockVector<double> &solution_delta
                                            ,const double delta_t)
{
    std::cout<<"Delta_t(inside newton):"<<delta_t<<std::endl;

    std::vector<BlockVector<double>> solution_delta_si(7,BlockVector<double>(dofs_per_block_m));
    std::vector<BlockVector<double>> current_sol_si(7,BlockVector<double>(dofs_per_block_m));
    
    current_sol_si[0]= solution_m;
    unsigned int new_iter_u = 0;
    
    //Staggered Iterations
    for(unsigned int k=1;k<8;++k)
    {
      get_error_residual(error_residual);

      if(k==1)
            error_residual_0 = error_residual;

      error_residual_norm = error_residual;
      error_residual_norm.normalize(error_residual_0);

      //checking
	    if( k>1 && (error_residual_norm.norm < param.pf.st_tol || error_residual.norm < param.pf.st_tol))
	    {      	
          //std::cout<<"Converged-Staggered Iterations"<<std::endl;
          solution_delta = solution_delta_si[k-1];
          break;
      }
    
      std::cout<<"Staggered Iteration:"<<k<<std::endl;
    

      ///////////////////////////////Linearising PhaseField-d//////////////////////////////////
      BlockVector<double> newton_update_d(dofs_per_block_m);

      error_residual.reset();
      error_residual_0.reset();
      error_residual_norm.reset();
      
      error_update.reset();
      error_update_0.reset();
      error_update_norm.reset();

      print_header_d();
      unsigned int new_iter_d = 0;

      for (; new_iter_d < param.newtonraphson.max_new_ite; ++new_iter_d)
      {

          std::cout << " " << std::setw(2) << new_iter_d << " " << std::flush;
    
          tangent_matrix_m = 0.0;
          system_rhs_m = 0.0;

          //Mesh induced initial crack + Phase field along the initial crack
          if(param.mod_strategy.strategy=="hybrid"){
          make_constraints_d(new_iter_d,param);
          }

          //Phase field induced initial crack
          if(param.mod_strategy.strategy=="phasefield"){        
          get_constrained_initial_d(new_iter_d,param);
          }

          current_sol_si[k-1] += newton_update_d; //Gives updated solution till now 
          assemble_system_d(param,current_sol_si[k-1],delta_t);	//Evaluation at solution calculated until now

          get_error_residual_d(error_residual);

          if(new_iter_d==0)
              error_residual_0 = error_residual;

          error_residual_norm = error_residual;
          error_residual_norm.normalize(error_residual_0);

          const std::pair<unsigned int,double>
                  lin_solver_output_d = solve_sys_d(param,newton_update_d);

          get_newton_update_error_d(newton_update_d,error_update);
          if(new_iter_d==0)
              error_update_0 = error_update;

          error_update_norm = error_update;
          error_update_norm.normalize(error_update_0);

          solution_delta_si[k] += newton_update_d; 

          //checking
          if(new_iter_d > 0 && (error_residual_norm.d < param.newtonraphson.res_tol_d
                                || error_residual.d < param.newtonraphson.res_tol_d))
          {      	
              std::cout<<"Converged-d"<<std::endl;
              print_footer_d();
              break;
          }

          std::cout << " | " << std::fixed << std::setprecision(3) << std::setw(4)
                              << std::scientific <<"        "<< lin_solver_output_d.first << "   "
                              << lin_solver_output_d.second 
                              << "  " << error_residual.d  << "  " << error_residual_norm.d
                              << "  " << error_update.d    << "  " << error_update_norm.d
                              << "  " << std::endl;

      }
      AssertThrow (new_iter_d < param.newtonraphson.max_new_ite,
                    ExcMessage("No convergence in nonlinear solver-d!"));
      
      ////////////////////////Linearising Displacement-u///////////////////////////////////////////
      BlockVector<double> newton_update_u(dofs_per_block_m);

      error_residual.reset();
      error_residual_0.reset();
      error_residual_norm.reset();
      
      error_update.reset();
      error_update_0.reset();
      error_update_norm.reset();

      print_header_u();
  //  unsigned int new_iter_u = 0;

      current_sol_si[k] = current_sol_si[k-1];

      for (; new_iter_u < param.newtonraphson.max_new_ite ; ++new_iter_u)
      {

          std::cout << " " << std::setw(2) << new_iter_u << " " << std::flush;
    
          tangent_matrix_m = 0.0;
          system_rhs_m = 0.0;

          make_constraints_u(new_iter_u,(delta_t/param.time.end_time),param);
          
          current_sol_si[k] += newton_update_u; //Gives updated solution till now
          
          assemble_system_u(param,current_sol_si[k]);	//Evaluation at solution calculated until now
  
          get_error_residual_u(error_residual);

          if(new_iter_u==0)
              error_residual_0 = error_residual;

          error_residual_norm = error_residual;
          error_residual_norm.normalize(error_residual_0);

          const std::pair<unsigned int,double>
                  lin_solver_output_u = solve_sys_u(newton_update_u);

          get_newton_update_error_u(newton_update_u,error_update);
          if(new_iter_u==0)
              error_update_0 = error_update;

          error_update_norm = error_update;
          error_update_norm.normalize(error_update_0);

          solution_delta_si[k] += newton_update_u;      
          //checking
          if(new_iter_u > 0 && (error_residual_norm.u < param.newtonraphson.res_tol_u || error_residual.u < param.newtonraphson.res_tol_u)
                          && (error_update_norm.u < param.newtonraphson.nu_tol_u || error_update.u < param.newtonraphson.nu_tol_u) )
          {      	
              std::cout<<"Converged-u"<<std::endl;
              print_footer_u();
              break;

          }

          const std::pair<double,double> energy = get_energy_p(param,newton_update_u);

          std::cout << " | " << std::fixed << std::setprecision(3) << std::setw(4)
                              << std::scientific << lin_solver_output_u.first << "   "
                              << lin_solver_output_u.second 
                              << "  " << error_residual.u  << "  " << error_residual_norm.u
                              << "  " << error_update.u    << "  " << error_update_norm.u
                              << "  " << energy.first      << "    " << energy.second
                              << "  " << std::endl;

      }
      if(new_iter_u==param.newtonraphson.max_new_ite)
      break;

    }
return new_iter_u;
}

template <int dim>
std::pair<unsigned int,double> Phasefield<dim>::solve_sys_d(const AllParameters &param,BlockVector<double> &newton_update)
{
  timer.enter_subsection("Linear solver-d");

  unsigned int lin_ite = 0;
  double lin_res = 0;

  SolverControl           solver_control (tangent_matrix_m.block(d_dof_m,d_dof_m).m(),
                                          param.linearsolver.cg_tol*system_rhs_m.block(d_dof_m).l2_norm());
  GrowingVectorMemory<Vector<double> > GVM;
  SolverCG<Vector<double>>    cg (solver_control,GVM);
  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(tangent_matrix_m.block(d_dof_m,d_dof_m), param.linearsolver.relax_prm);
  cg.solve (tangent_matrix_m.block(d_dof_m,d_dof_m), newton_update.block(d_dof_m), system_rhs_m.block(d_dof_m),preconditioner);


  lin_ite = solver_control.last_step();
  lin_res = solver_control.last_value();

  constraints_m.distribute (newton_update);		//?? //no effect
  timer.leave_subsection();

 return std::make_pair(lin_ite,lin_res);
}
  

template <int dim>
std::pair<unsigned int,double> Phasefield<dim>::solve_sys_u(BlockVector<double> &newton_update)
{
  timer.enter_subsection("Linear solver-u");
  
  SparseDirectUMFPACK k_uu;
  k_uu.initialize(tangent_matrix_m.block(u_dof_m,u_dof_m)); 
  k_uu.vmult(newton_update.block(u_dof_m),system_rhs_m.block(u_dof_m)); 

  unsigned int lin_ite = 0;
  double lin_res = 0;

  constraints_m.distribute (newton_update);

  timer.leave_subsection();
 
  return std::make_pair(lin_ite,lin_res);
}


template <int dim>
void Phasefield<dim>::make_constraints_u(unsigned int &itr,const double load_ratio,const AllParameters &param){
    std::cout<<" CST "<<std::flush;

    if(itr>1){				//no effect of removing
        return;
    }

   constraints_m.clear();

   DoFTools::make_hanging_node_constraints (dof_handler_m,
                                           constraints_m);

   std::vector<bool> component_mask(dim+1, false);
     
   //Tension test
   if(param.test_case.test == "tension")
   {
     if(param.mod_strategy.comp_strategy=="normal")
     {
      { 
      std::vector<bool> component_mask1(dim+1, false);

      if(param.bc.uxb == "fixed")
        component_mask1[0] = true;
      else if(param.bc.uxb == "free")
        component_mask1[0] = false;
      
      component_mask1[1] = true;
      component_mask1[2] = false;

      VectorTools::interpolate_boundary_values(dof_handler_m, 3,
                                                  ZeroFunction<dim>(n_components_m), constraints_m, component_mask1);
      }      
      {
      std::vector<bool> component_mask2(dim+1, false);

      if(param.bc.uxt == "fixed")
        component_mask2[0] = true;
      else if(param.bc.uxt == "free")
        component_mask2[0] = false;
      
      component_mask2[1] = true;
      component_mask2[2] = false;

      VectorTools::interpolate_boundary_values(dof_handler_m, 4,
                                  BoundaryTension<dim>(itr,load_ratio,param.pf.u_total), constraints_m, component_mask2);
      }
     }
      
     if(param.mod_strategy.comp_strategy=="lefm")
      {         
        
        std::vector<bool> component_mask(dim+1, false);

        component_mask[0] = true;
        component_mask[1] = true;
        component_mask[2] = false;

        VectorTools::interpolate_boundary_values(dof_handler_m, 1,
                                    Reference_solution<dim>(itr,param.mod_strategy.steps_ft
                                                            ,param.pf.g_c
                                                            ,param.materialmodel.lambda
                                                            ,param.materialmodel.mu), constraints_m, component_mask);
        VectorTools::interpolate_boundary_values(dof_handler_m, 2,
                                    Reference_solution<dim>(itr,param.mod_strategy.steps_ft
                                                            ,param.pf.g_c
                                                            ,param.materialmodel.lambda
                                                            ,param.materialmodel.mu), constraints_m, component_mask);
        VectorTools::interpolate_boundary_values(dof_handler_m, 3,
                                    Reference_solution<dim>(itr,param.mod_strategy.steps_ft
                                                            ,param.pf.g_c
                                                            ,param.materialmodel.lambda
                                                            ,param.materialmodel.mu), constraints_m, component_mask);
        VectorTools::interpolate_boundary_values(dof_handler_m, 4,
                                    Reference_solution<dim>(itr,param.mod_strategy.steps_ft
                                                            ,param.pf.g_c
                                                            ,param.materialmodel.lambda
                                                            ,param.materialmodel.mu), constraints_m, component_mask);                                                                                                                                                            
                                                            
        
      }
      
              
   }

//Shear Test
   else if(param.test_case.test == "shear")
   {
      { 
      component_mask[0] = true;
      component_mask[1] = true;
      VectorTools::interpolate_boundary_values(dof_handler_m, 3,
                                                  ZeroFunction<dim>(n_components_m), constraints_m, component_mask);
      }      
      {
      component_mask[0] = true;
      component_mask[1] = true;
      VectorTools::interpolate_boundary_values(dof_handler_m, 4,
                                                  BoundaryShear<dim>(itr,load_ratio,param.pf.u_total), constraints_m, component_mask);
      }
      { 
      component_mask[0] = false;
      component_mask[1] = true;
      VectorTools::interpolate_boundary_values(dof_handler_m, 1,
                                                  ZeroFunction<dim>(n_components_m), constraints_m, component_mask);
      }      
      {
      component_mask[0] = false;
      component_mask[1] = true;
      VectorTools::interpolate_boundary_values(dof_handler_m, 2,
                                                  ZeroFunction<dim>(n_components_m), constraints_m, component_mask);
      }

   }
   
   constraints_m.close ();
}

template <int dim>
void Phasefield<dim>::make_constraints_d(unsigned int &itr,const AllParameters &param){

    if(itr>1){				//no effect of removing
        return;
    }

   constraints_m.clear();

   DoFTools::make_hanging_node_constraints (dof_handler_m,
                                            constraints_m);
   //Tension test
   if(param.test_case.test == "tension")
   {
      //d=1 on slit DBC 
      if(current_time_m==param.time.delta_t){
      std::vector<bool> component_mask(dim+1, false);

      component_mask[0] = false;
      component_mask[1] = false;
      component_mask[2] = true;
      VectorTools::interpolate_boundary_values(dof_handler_m, 0,
                                                      InitialCrack<dim>(itr), constraints_m, component_mask);
      }
      else
      {
      std::vector<bool> component_mask(dim+1, false);

      component_mask[0] = false;
      component_mask[1] = false;
      component_mask[2] = true;
      VectorTools::interpolate_boundary_values(dof_handler_m, 0,
                                                  ZeroFunction<dim>(n_components_m), constraints_m, component_mask);
      }
      
       
   }
   
//Shear Test
   else if(param.test_case.test == "shear")
   {
      
      // {
      // std::vector<bool> component_mask(dim+1, false);

      // component_mask[0] = false;
      // component_mask[1] = false;
      // component_mask[2] = true;
      // VectorTools::interpolate_boundary_values(dof_handler_m, 0,
      //                                                 InitialCrack<dim>(itr), constraints_m, component_mask);
      // }
   }
   
   //constraints_m.close ();
}

template <int dim>
void Phasefield<dim>::print_header_d(){
    const unsigned int l_width = 80;
        for (unsigned int i = 0; i < l_width; ++i)
        {
            std::cout << "_";
        }
        std::cout << std::endl;

        std::cout << "SOLVER STEP"
                  << "| LIN_IT  LIN_RES     "
		              << "RES_D     RES_D_N     "                    
		              << "NU_D      NU_D_N      "
		              << std::endl;

        for (unsigned int i = 0; i < l_width; ++i)
        {
            std::cout << "_";
        }
        std::cout << std::endl;

}
template <int dim>
void Phasefield<dim>::print_header_u(){
    const unsigned int l_width = 105;
        for (unsigned int i = 0; i < l_width; ++i)
        {
            std::cout << "_";
        }
        std::cout << std::endl;

        std::cout << "SOLVER STEP"
                  << "| LIN_IT  LIN_RES     "
                  << "RES_U    RES_U_N    "                    
                  << "NU_U       NU_U_N   "
                  << "Total Energy  Fracture Energy"
                  << std::endl;

        for (unsigned int i = 0; i < l_width; ++i)
        {
            std::cout << "_";
        }
        std::cout << std::endl;

}
template <int dim>
void Phasefield<dim>::print_footer_d(){
    const unsigned int l_width = 80;
        for (unsigned int i = 0; i < l_width; ++i)
        {
            std::cout << "_";
        }
        std::cout << std::endl;

        std::cout << "Errors:" << std::endl
                  << "Res_d_n: \t\t" << error_residual_norm.d << std::endl
		  << std::endl;

}

template <int dim>
void Phasefield<dim>::print_footer_u(){
    const unsigned int l_width = 105;
        for (unsigned int i = 0; i < l_width; ++i)
        {
            std::cout << "_";
        }
        std::cout << std::endl;

        std::cout << "Errors:" << std::endl
                  << "Res_u_n: \t\t" << error_residual_norm.u << std::endl
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
}

template <int dim>
void Phasefield<dim>::get_error_residual_d(Error& error_residual){

    BlockVector<double> err_res(dofs_per_block_m);
    for (unsigned int i=0;i<dof_handler_m.n_dofs();++i) {
        if(!constraints_m.is_constrained(i)){
            err_res(i)=system_rhs_m(i);
           }

        }
    error_residual.norm = err_res.l2_norm();
    error_residual.d = err_res.block(d_dof_m).l2_norm();
}

template <int dim>
void Phasefield<dim>::get_error_residual_u(Error& error_residual){

    BlockVector<double> err_res(dofs_per_block_m);
    
    for (unsigned int i=0;i<dof_handler_m.n_dofs();++i) {
        if(!constraints_m.is_constrained(i)){
            err_res(i)=system_rhs_m(i);
           }
        }
    error_residual.norm = err_res.l2_norm();
    error_residual.u = err_res.block(u_dof_m).l2_norm();
}

template <int dim>
void Phasefield<dim>::get_newton_update_error_d(const BlockVector<double> &newton_update
					                                      ,Error& error_update){

    BlockVector<double> err_update(dofs_per_block_m);
    for (unsigned int i=0;i<dof_handler_m.n_dofs();++i) {
        if(!constraints_m.is_constrained(i)){
            err_update(i)= newton_update(i);
           }
        }
    error_update.norm = err_update.l2_norm();
    error_update.d = err_update.block(d_dof_m).l2_norm();
}

template <int dim>
void Phasefield<dim>::get_newton_update_error_u(const BlockVector<double> &newton_update
					     ,Error& error_update){

    BlockVector<double> err_update(dofs_per_block_m);
    for (unsigned int i=0;i<dof_handler_m.n_dofs();++i) {
        if(!constraints_m.is_constrained(i)){
            err_update(i)= newton_update(i);
           }
        }
    error_update.norm = err_update.l2_norm();
    error_update.u = err_update.block(u_dof_m).l2_norm();
}

template <int dim>
void Phasefield<dim>::import_mesh(const AllParameters &param){

    std::string grid_name;
    grid_name += param.geometrymodel.meshfile;

    GridIn<dim> grid_in;
    grid_in.attach_triangulation(triangulation_m);
    std::ifstream input_file(grid_name.c_str());
    Assert(dim==2, ExcInternalError());
    grid_in.read_abaqus(input_file,false);
//    grid_in.read_ucd(input_file,false);

    GridTools::scale(param.geometrymodel.grid_scale,triangulation_m);
    triangulation_m.refine_global (param.geometrymodel.gl_ref);

    const bool write_grid = true;
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

    set_boundary_id(param);
}

template <int dim>
void Phasefield<dim>::output_results (const AllParameters &param,unsigned int cycle) const
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
    
//    MappingQEulerian<dim> q_mapping(param.fesys.fe_degree, dof_handler_m, soln);
      MappingQ <dim> q_mapping(param.fesys.fe_degree);

    data_out.build_patches(q_mapping, param.fesys.fe_degree);
    std::string filename ("solution-"
                          + std::to_string(dim)
                          + "d-"
                          + std::to_string(cycle/param.time.op_freq)
                          + ".vtu");
    std::ofstream output(filename.c_str());
    data_out.write_vtu(output);

    std::vector<std::pair<double,std::string>>  times_and_names;
    times_and_names.push_back(std::make_pair(current_time_m,filename));
    std::ofstream pvd_output("solution.pvd");
    DataOutBase::write_pvd_record(pvd_output,times_and_names);
  
}

  

template <int dim>
void Phasefield<dim>::run(const AllParameters &param){


    using namespace constants;

    BlockVector<double>       solution_delta(dofs_per_block_m);

    long double present_time_tol;
    double delta_t;

  //  if(param.mod_strategy.comp_strategy=="normal")
    delta_t = param.time.delta_t;
    
  //  else if(param.mod_strategy.comp_strategy=="lefm")
  //  delta_t = (1/(param.mod_strategy.steps_ft*param.mod_strategy.fac_ft));

    present_time_tol = time_tol_fac * delta_t;
    current_time_m = param.time.start_time;
    unsigned int current_timestep = 0;

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
                if(param.mod_strategy.comp_strategy=="normal")
                {
                  if (cell_vertex[0] <= 1.0 && cell_vertex[0] >= 0 /* 0.48*/
                      && cell_vertex[1] <= 0.505 && cell_vertex[1] >= 0.495)
                  {
                      cell->set_refine_flag();
                      break;
                  }
                }
                else if(param.mod_strategy.comp_strategy=="lefm")
                {
                  if (cell_vertex[0] <= 0.5 && cell_vertex[0] >= -0.1 
                      && cell_vertex[1] <= 0.015 && cell_vertex[1] >= -0.015)
                  {
                      cell->set_refine_flag();
                      break;
                  }
                }
	            }
          }

    triangulation_m.execute_coarsening_and_refinement();
    setup_system();
    setup_qph();
    solution_delta.reinit(dofs_per_block_m);
    }

    //Phase-field Induced crack
    if(param.mod_strategy.strategy=="phasefield")
    {
      double min_cell_diameter = GridTools::minimal_cell_diameter(triangulation_m);
      extract_initialcrack_d_index(min_cell_diameter,param);
    }

    output_results(param,current_timestep);
    statistics.add_value("Time",current_time_m);
    statistics.set_precision("Time",6);
    compute_load(param,solution_m);
    get_energy_v(param,solution_m);

    if(param.mod_strategy.comp_strategy=="lefm")
    compute_KI(param,current_timestep);

  //  if(param.mod_strategy.comp_strategy=="normal")
    current_time_m += param.time.delta_t;
    
  //  else if(param.mod_strategy.comp_strategy=="lefm")
  //  current_time_m += delta_t;

    current_timestep++;

/////////Printing////////
	std::cout<<"Min cell dia:"<<GridTools::minimal_cell_diameter(triangulation_m)<<std::endl;
	std::cout<<"Max cell dia:"<<GridTools::maximal_cell_diameter(triangulation_m)<<std::endl;
	std::cout<<"Parameter Global Ref: "<<param.geometrymodel.gl_ref<<std::endl;
  std::cout<<"Parameter Locl Ref: "<<param.geometrymodel.lc_ref<<std::endl;
  std::cout<<"Parameter Grid scale: "<<param.geometrymodel.grid_scale<<std::endl;
  std::cout<<"Parameter lambda: "<<param.materialmodel.lambda<<std::endl;
	std::cout<<"Parameter mu: "<<param.materialmodel.mu<<std::endl;
	std::cout<<"Parameter vis: "<<param.materialmodel.viscosity<<std::endl;
	std::cout<<"Parameter g_c: "<<param.pf.g_c<<std::endl;
	std::cout<<"Parameter l: "<<param.pf.l<<std::endl;
	std::cout<<"Parameter k: "<<param.pf.k<<std::endl;
	std::cout<<"Parameter delta_t: "<</*param.time.*/delta_t<<std::endl;
	std::cout<<"Parameter end_time: "<<param.time.end_time<<std::endl;
	std::cout<<"Parameter u_total: "<<param.pf.u_total<<std::endl;
	std::cout<<"Parameter time change interval: "<<param.time.time_change_interval<<std::endl;
  std::cout<<"Parameter Newton: max iterations: "<<param.newtonraphson.max_new_ite<<std::endl;
	std::cout<<"Parameter Newton: res tol_u: "<<param.newtonraphson.res_tol_u<<std::endl;
	std::cout<<"Parameter Newton: res tol_d: "<<param.newtonraphson.res_tol_d<<std::endl;
	std::cout<<"Parameter Newton: new_upd tol_u: "<<param.newtonraphson.nu_tol_u<<std::endl;
	std::cout<<"Parameter Newton: new_upd tol_d: "<<param.newtonraphson.nu_tol_d<<std::endl;
	std::cout<<std::endl;

///////////////////////
     

    unsigned int new_iter=0;
    double prev_delta_t = 0; 

    while (current_time_m < param.time.end_time + present_time_tol)
	 
    {		
      std::cout<<std::endl;
      std::cout<<std::setprecision(5)<<"Current_time:"<<current_time_m<<std::endl;
      std::cout<<"Current time step: "<<current_timestep<<std::endl;

      solution_delta = 0.0;
      new_iter = solve_nonlinear_newton(param,solution_delta,delta_t);

      if(new_iter == param.newtonraphson.max_new_ite)
      {
          prev_delta_t = delta_t;
          delta_t = param.time.alpha*delta_t;
          if(delta_t < param.time.beta*param.time.delta_t)
          {
            std::cout<<"delta_t:"<<delta_t<<std::endl;
            std::cout<<"param.time.delta_t:"<<param.time.delta_t<<std::endl;
            std::cout<<"beta*param.time.delta_t:"<<param.time.beta*param.time.delta_t<<std::endl;
            std::cout<<std::endl;
            std::cout << "No convergence in prescribed limits!" << std::endl;
            std::cerr << "No convergence in prescribed limits!" << std::endl
                      << "Aborting!" << std::endl;
                      exit(0);
          }
          std::cout<<"Solver did not converge! Adjusting time steps to "<<delta_t<<std::endl;
          current_time_m -= prev_delta_t;
          current_time_m += delta_t;
      }
      else
      {

        solution_m += solution_delta;	
        old_solution_m = solution_m;

        statistics.add_value("Time",current_time_m);
        statistics.set_precision("Time",6);
        compute_load(param,solution_m);
        get_energy_v(param,solution_m);
        if(param.mod_strategy.comp_strategy=="lefm")
        compute_KI(param,current_timestep);

        if((current_timestep%param.time.op_freq) == 0 ||
            current_time_m==param.time.delta_t ||/*
            std::fabs(current_time_m - 0.0046)<1e-8 ||
            std::fabs(current_time_m - 0.00465)<1e-8 ||
            std::fabs(current_time_m - 0.0047)<1e-8 ||
            std::fabs(current_time_m - 0.00475)<1e-8 ||
            std::fabs(current_time_m - 0.0048)<1e-8 ||
            std::fabs(current_time_m - 0.00485)<1e-8 ||
            std::fabs(current_time_m - 0.0049)<1e-8 ||
            std::fabs(current_time_m - 0.00495)<1e-8 ||
            std::fabs(current_time_m - 0.0050)<1e-8 ||
            std::fabs(current_time_m - 0.00505)<1e-8 ||
            std::fabs(current_time_m - 0.0051)<1e-8 ||
            std::fabs(current_time_m - 0.00515)<1e-8 ||
            std::fabs(current_time_m - 0.0052)<1e-8 ||
            std::fabs(current_time_m - 0.00525)<1e-8 ||
            std::fabs(current_time_m - 0.0053)<1e-8 ||
            std::fabs(current_time_m - 0.00535)<1e-8 ||
            std::fabs(current_time_m - 0.0054)<1e-8 ||
            std::fabs(current_time_m - 0.00545)<1e-8 ||
            std::fabs(current_time_m - 0.0055)<1e-8 ||
            std::fabs(current_time_m - 0.00555)<1e-8 ||
            std::fabs(current_time_m - 0.0056)<1e-8 ||
            std::fabs(current_time_m - 0.00565)<1e-8 ||
            std::fabs(current_time_m - 0.0057)<1e-8 ||
            std::fabs(current_time_m - 0.00575)<1e-8 ||
            std::fabs(current_time_m - 0.0058)<1e-8 ||
            std::fabs(current_time_m - 0.00585)<1e-8 ||
            std::fabs(current_time_m - 0.0059)<1e-8 ||
            std::fabs(current_time_m - 0.00595)<1e-8 ||
            std::fabs(current_time_m - 0.0060)<1e-8 ||
            std::fabs(current_time_m - 0.00605)<1e-8 ||
            std::fabs(current_time_m - 0.0061)<1e-8 ||
            std::fabs(current_time_m - 0.0062)<1e-8 ||
            std::fabs(current_time_m - 0.0063)<1e-8 ||
            std::fabs(current_time_m - 0.0080)<1e-8 ||
            std::fabs(current_time_m - 0.0085)<1e-8 ||
            std::fabs(current_time_m - 0.0090)<1e-8 ||
            std::fabs(current_time_m - 0.0095)<1e-8 ||
            std::fabs(current_time_m - 0.0100)<1e-8 ||
            std::fabs(current_time_m - 0.0105)<1e-8 ||
            std::fabs(current_time_m - 0.0110)<1e-8 ||
            std::fabs(current_time_m - 0.0115)<1e-8 ||
            std::fabs(current_time_m - 0.0120)<1e-8 ||
            std::fabs(current_time_m - 0.0125)<1e-8 ||
            std::fabs(current_time_m - 0.0130)<1e-8 ||
            std::fabs(current_time_m - 0.0134)<1e-8 ||
            std::fabs(current_time_m - 0.0135)<1e-8 ||
            std::fabs(current_time_m - 0.0140)<1e-8 ||
            std::fabs(current_time_m - 0.0145)<1e-8 ||
            std::fabs(current_time_m - 0.0148)<1e-8 ||
            std::fabs(current_time_m - 0.0150)<1e-8 ||
            std::fabs(current_time_m - 0.0155)<1e-8*/
            compute_end_load(param,solution_m) < 1.0
             )//Tension
        {
          std::cout<<"time:"<<current_time_m<<std::endl;
              output_results(param,current_timestep);
        }
     
        std::ofstream stat_file ("statistics_m");
        statistics.write_text (stat_file,
                                    TableHandler::simple_table_with_separate_column_description);
        stat_file.close();
        if(param.test_case.test == "tension" && param.mod_strategy.comp_strategy=="normal" )
        {
          if(current_time_m >= param.time.time_change_interval)
          {
            static bool once=false;
            if(!once)
            {
              delta_t = param.time.delta_t_f;
              once =true;
            }
          }
        }
        if(param.mod_strategy.comp_strategy=="lefm")
        {
          if(current_timestep>param.mod_strategy.fac_ft*param.mod_strategy.steps_ft)
          {
            std::cout<<"Reached end of the timesteps"<<std::endl;
            break;
          }
        }
        current_time_m += delta_t;
        ++current_timestep;
      }
               
   }

}

template class thesis::Phasefield<2>;
template class thesis::Phasefield<3>;

