#include "../include/Phasefield.h"
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
{
    import_mesh(param);
    dof_handler_m.initialize(triangulation_m, fe_m);
    setup_system();
    determine_comp_extractor();
}



template <int dim>
Phasefield<dim>::~Phasefield ()
{
  dof_handler_m.clear ();
}


template <int dim>
void Phasefield<dim>::setup_system ()
{

  dof_handler_m.distribute_dofs (fe_m);
  DoFRenumbering::block_wise(dof_handler_m); 
  DoFTools::count_dofs_per_block (dof_handler_m, dofs_per_block_m);

  std::cout << "   Number of active cells:       "
                            << triangulation_m.n_active_cells()
                            << std::endl;
  std::cout << "   Number of degrees of freedom: "
                            << dof_handler_m.n_dofs()
                            << std::endl;

  constraints_m.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler_m,
                                           constraints_m);
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
std::cout<<"setup solution.size: "<<solution_m.size()<<std::endl;

  solution_m.reinit (dofs_per_block_m);
  solution_m.collect_sizes();
std::cout<<"setup solution.size: "<<solution_m.size()<<std::endl;
std::cout<<"setup solution_u.size: "<<solution_m.block(u_dof_m).size()<<std::endl;
std::cout<<"setup solution_d.size: "<<solution_m.block(d_dof_m).size()<<std::endl;
for(auto &c:dofs_per_block_m)
	std::cout<<c<<"\t";
 old_solution_m.reinit (dofs_per_block_m);
 old_solution_m.collect_sizes();
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
	      }
          }

}


template <int dim>
void Phasefield<dim>::assemble_system (const AllParameters &param,BlockVector<double> & newton_update)

{

  FEValues<dim> fe_values (fe_m, quadrature_formula_m,
                           update_values   | update_gradients |
                           update_quadrature_points | update_JxW_values);
  const unsigned int   dofs_per_cell = fe_m.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula_m.size();
  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  double History;

  for (const auto &cell : dof_handler_m.active_cell_iterators())
    {
      cell_matrix = 0;

      fe_values.reinit(cell);

      cell->get_dof_indices (local_dof_indices);

      std::vector<SymmetricTensor<2,dim>> epsilon_vals(n_q_points);
      fe_values[u_extractor].get_function_symmetric_gradients(newton_update,epsilon_vals);
      
///Testing///
      BlockVector<double> tmp;
      tmp.reinit(dofs_per_block_m);
      tmp.collect_sizes();
      std::vector<double> d_vals(n_q_points);
      fe_values[d_extractor].get_function_values(/*newton_update*/tmp,d_vals);
      
        for (unsigned int q = 0; q < n_q_points; ++q){
		
	    if(current_time_m == param.time.delta_t)
		History = 0;
	    else
		History = old_history_m;
	    const SymmetricTensor<2,dim> sigma_plus = get_stress_plus(param.materialmodel.lambda,param.materialmodel.mu,epsilon_vals[q]);
	    const SymmetricTensor<4,dim> BigC_plus = get_BigC_plus(param.materialmodel.lambda,param.materialmodel.mu,epsilon_vals[q]);
	    const SymmetricTensor<4,dim> BigC_minus = get_BigC_minus(param.materialmodel.lambda,param.materialmodel.mu,epsilon_vals[q]);
//Printing/////////////////
/*	    if(cell==dof_handler_m.begin_active() && q==2){
		std::cout<<std::endl;
		std::cout<<"d_vals[q] _system "<<d_vals[q]<<std::endl;
//		std::cout<<"newton update norm_system : "<<newton_update.l2_norm()<<std::endl;
//		}

		std::cout<<"History_assemble_system: "<<History<<std::endl;
		std::cout<<"sigma_plus_assem_sys:"<<std::endl;
		print_tensor(sigma_plus);
		std::cout<<"BigC_plus:"<<std::endl;
		print_tensor(BigC_plus);
		std::cout<<"BigC_minus:"<<std::endl;
		print_tensor(BigC_minus);

		}
*/
//////////////////////////

            for (unsigned int i = 0; i < dofs_per_cell; ++i){

                const Tensor<2, dim> grad_shape_i_u = fe_values[u_extractor].gradient(i, q);
                const Tensor<1, dim> grad_shape_i_d = fe_values[d_extractor].gradient(i, q);
		const double shape_i_d = fe_values[d_extractor].value(i,q);


                for (unsigned int j = 0; j < dofs_per_cell; ++j){

                const SymmetricTensor<2, dim> sym_grad_shape_j_u = fe_values[u_extractor].symmetric_gradient(j, q);
                const Tensor<1, dim> grad_shape_j_d = fe_values[d_extractor].gradient(j, q);
		const double shape_j_d = fe_values[d_extractor].value(j,q);

			if((dof_block_identifier_m[i] == d_dof_m) && (dof_block_identifier_m[j] == d_dof_m))
			{
                 	 cell_matrix(i, j) += ( (param.pf.g_c/param.pf.l)*shape_i_d*shape_j_d
				              + (param.pf.g_c*param.pf.l)*scalar_product(grad_shape_i_d,grad_shape_j_d)
					      + 2*shape_i_d*shape_j_d*History
					      + (param.materialmodel.viscosity/param.time.delta_t)*shape_i_d*shape_j_d ) *fe_values.JxW(q);
			}
			else if((dof_block_identifier_m[i] == d_dof_m) && (dof_block_identifier_m[j] == u_dof_m))
			{
			cell_matrix(i,j) += 0;
			}
			else if((dof_block_identifier_m[i] == u_dof_m) && (dof_block_identifier_m[j] == d_dof_m))
			{
			cell_matrix(i,j) += ( 2*scalar_product(grad_shape_i_u,sigma_plus)*(1-d_vals[q])*shape_j_d )*fe_values.JxW(q);
			}
			else if((dof_block_identifier_m[i] == u_dof_m) && (dof_block_identifier_m[j] == u_dof_m))
			{
			SymmetricTensor<2,dim> tmp1,tmp2;
			double_contract(tmp1, BigC_plus, sym_grad_shape_j_u);
			double_contract(tmp2, BigC_minus,sym_grad_shape_j_u);

			cell_matrix(i,j) += (- scalar_product(grad_shape_i_u,tmp1)*(1 - d_vals[q])*(1 - d_vals[q])
					     - param.pf.k*scalar_product(grad_shape_i_u,tmp1)
					     - scalar_product(grad_shape_i_u,tmp2)  )*fe_values.JxW(q);
			}
                }
            }
          }
/*
        for(unsigned int i=0;i<dofs_per_cell;++i)
            for(unsigned int j=0;j<i;++j)
                cell_matrix(i,j)=cell_matrix(j,i);
*/  
//     constraints_m.distribute_local_to_global(cell_matrix,cell_rhs,
//                                  local_dof_indices,tangent_matrix_m,system_rhs_m,false);  //copy local to global
     constraints_m.distribute_local_to_global(cell_matrix,local_dof_indices,tangent_matrix_m);  //copy local to global


  }
}


template <int dim>
void Phasefield<dim>::assemble_rhs(const AllParameters &param,BlockVector<double> & newton_update)

{

  FEValues<dim> fe_values (fe_m, quadrature_formula_m,
                           update_values   | update_gradients |
                           update_quadrature_points | update_JxW_values);
  
  FEFaceValues<dim> fe_values_face(fe_m, face_quadrature_formula_m,
                                   update_values | update_quadrature_points |
                                   update_normal_vectors | update_JxW_values);

  const unsigned int   dofs_per_cell = fe_m.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula_m.size();
  const unsigned int   n_face_q_points    = face_quadrature_formula_m.size(); 
  
  Vector<double>       cell_rhs (dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  std::vector<Vector<double> >      traction_force_values (n_face_q_points,Vector<double>(dim+1));
  const BoundaryForce<dim> 		    traction_force;

  double History;
  Tensor<1,dim> load_value;
  Tensor<1,dim> load_value_1;
   
    
  for (const auto &cell : dof_handler_m.active_cell_iterators())
    {

	cell_rhs = 0;

	fe_values.reinit(cell);
	cell->get_dof_indices (local_dof_indices);

	Tensor<1,dim> body_force;
	body_force = 0;
	
	std::vector<SymmetricTensor<2,dim>> epsilon_vals(n_q_points);
	fe_values[u_extractor].get_function_symmetric_gradients(newton_update,epsilon_vals);
 
//Testing///
        BlockVector<double> tmp;
        tmp.reinit(dofs_per_block_m);
        tmp.collect_sizes();
 
	std::vector<Tensor<1,dim>> grad_d(n_q_points);
        fe_values[d_extractor].get_function_gradients(/*newton_update*/tmp,grad_d);

 	std::vector<double> d_vals(n_q_points);
	fe_values[d_extractor].get_function_values(/*newton_update*/tmp,d_vals);
  
	std::vector<double> old_d_vals(n_q_points);
	fe_values[d_extractor].get_function_values(/*old_solution_m*/tmp,old_d_vals);
///Testing/////  
	std::vector<Tensor<1,dim>> old_u_vals(n_q_points);
	fe_values[u_extractor].get_function_values(old_solution_m,old_u_vals);
/////////////// 
        for (unsigned int q = 0; q < n_q_points; ++q){
        		
	    if(current_time_m == param.time.delta_t)
		History = 0;
	    else
		History = old_history_m;
	  
		const SymmetricTensor<2,dim> sigma_plus = get_stress_plus(param.materialmodel.lambda,param.materialmodel.mu,epsilon_vals[q]);
		const SymmetricTensor<2,dim> sigma_minus= get_stress_minus(param.materialmodel.lambda,param.materialmodel.mu,epsilon_vals[q]);
//Printing////////
//		std::cout<<"cell"<<cell<<std::endl;
/*	    if(cell==dof_handler_m.begin_active() && q==0){
		std::cout<<std::endl;
		std::cout<<"d_vals[q] _rhs "<<d_vals[q]<<std::endl;
//		std::cout<<"newton update norm_rhs : "<<newton_update.l2_norm()<<std::endl;
		std::cout<<"old_d_vals[q] _rhs : "<<old_d_vals[q]<<std::endl;
		std::cout<<"old_u_vals[q] _rhs : "<<old_u_vals[q]<<std::endl;
		std::cout<<"History_assem_rhs: "<<History<<std::endl;
//		std::cout<<"sigma_plus_assem_rhs:"<<std::endl;
		print_tensor(sigma_plus);
		std::cout<<"sigma_minus_assem_rhs:"<<std::endl;
		print_tensor(sigma_minus);
		std::cout<<"grad_d[q]:"<<grad_d[q]<<std::endl;
		}
*/
//////////////////////////

      
		  for (unsigned int i = 0; i < dofs_per_cell; ++i){
		   
		  const Tensor<2, dim> grad_shape_i_u = fe_values[u_extractor].gradient(i, q);
		  const Tensor<1, dim> shape_i_u = fe_values[u_extractor].value(i, q);
		  const Tensor<1, dim> grad_shape_i_d = fe_values[d_extractor].gradient(i, q);
		  const double shape_i_d = fe_values[d_extractor].value(i,q);
	
			if(dof_block_identifier_m[i] == d_dof_m )
                        {
                         cell_rhs(i) -= ( (param.pf.g_c/param.pf.l)*shape_i_d *d_vals[q]
					+ (param.pf.g_c*param.pf.l)*scalar_product(grad_shape_i_d,grad_d[q])
					- 2*shape_i_d*(1 - d_vals[q]) *History
					+ (param.materialmodel.viscosity/param.time.delta_t)*shape_i_d * d_vals[q]
					- (param.materialmodel.viscosity/param.time.delta_t)*shape_i_d * old_d_vals[q] 
					   )*fe_values.JxW(q);

                        }
                        else if(dof_block_identifier_m[i] == u_dof_m)
                        {
                        cell_rhs(i) -= (-scalar_product(grad_shape_i_u,sigma_plus)*(1 - d_vals[q])*(1-d_vals[q])
					- param.pf.k*scalar_product(grad_shape_i_u,sigma_plus)
					- scalar_product(grad_shape_i_u,sigma_minus)
					+ scalar_product(shape_i_u,body_force)
					 ) *fe_values.JxW(q);
                        }

           	}
	}

	
	for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
	{
	          if (cell->face(face)->at_boundary() && cell->face(face)->boundary_id() == 3)
            	  {
              		fe_values_face.reinit(cell, face);
              		traction_force.vector_value_list(fe_values_face.get_quadrature_points(),
                                               traction_force_values);
              		
			for (unsigned int q_point = 0; q_point < n_face_q_points;++q_point)
                	{
//TO -traction                  	
				Tensor<1, dim> rhs_values;
				rhs_values = 0;

				for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    		cell_rhs(i) += (fe_values_face[u_extractor].value(i, q_point) * rhs_values
                                	       * fe_values_face.JxW(q_point));
				//For Printing Load values		
				Tensor<2,dim> stress_display;
				double tr_eps;
				tr_eps = trace(epsilon_vals[q_point]);
				stress_display = param.materialmodel.lambda*tr_eps*unit_symmetric_tensor<dim>()
						+ 2*param.materialmodel.mu*epsilon_vals[q_point];
				load_value += stress_display*fe_values_face.normal_vector(q_point)*fe_values_face.JxW(q_point);
				///////////////////////////
			}
            	  }
		/*
		if(cell->face(face)->at_boundary() && (cell->face(face)->boundary_id() == 0 ||
							cell->face(face)->boundary_id() == 1 ||
							cell->face(face)->boundary_id() == 2) )
		{
    	       		fe_values_face.reinit(cell, face);
				for (unsigned int q_point = 0; q_point < n_face_q_points;++q_point)
                	{

				Tensor<2,dim> stress_display_1;
				double tr_eps;
				tr_eps = trace(epsilon_vals[q_point]);
				stress_display_1 = param.materialmodel.lambda*tr_eps*unit_symmetric_tensor<dim>()
						+ 2*param.materialmodel.mu*epsilon_vals[q_point];
				load_value_1 += stress_display_1*fe_values_face.normal_vector(q_point)*fe_values_face.JxW(q_point);
			}			
		}*/
	}
//	cell->get_dof_indices (local_dof_indices);
    	constraints_m.distribute_local_to_global(cell_rhs,
                                  local_dof_indices,
                                  system_rhs_m);//copy local to global
  
	constraints_m.set_zero(system_rhs_m);
  }

std::cout<<std::endl;
std::cout<<"Boundary3_Load_x: "<<load_value[0]<<std::endl;
std::cout<<"Boundary3_Load_y: "<<load_value[1]<<std::endl;
//std::cout<<"Boundary012_Load_x: "<<load_value_1[0]<<std::endl;
//std::cout<<"Boundary012_Load_y: "<<load_value_1[1]<<std::endl;

}

template <int dim>
void Phasefield<dim>::solve_nonlinear_newton(const AllParameters &param,
                                                 BlockVector<double> &solution_delta){

    BlockVector<double> newton_update(dofs_per_block_m);

    error_residual.reset();
    error_residual_0.reset();
    error_residual_norm.reset();

    print_header();
    unsigned int new_iter = 0;
    for (; new_iter < param.newtonraphson.max_new_ite; ++new_iter) {

        std::cout << " " << std::setw(2) << new_iter << " " << std::flush;

        tangent_matrix_m = 0.0;
        system_rhs_m = 0.0;

        make_constraints(new_iter,current_time_m);

        assemble_system(param,newton_update/*solution_delta*/);
        assemble_rhs(param,newton_update/*solution_delta*/);

        get_error_residual(error_residual);

        if(new_iter==0)
            error_residual_0 = error_residual;

        error_residual_norm = error_residual;
        error_residual_norm.normalize(error_residual_0);

        if(new_iter > 0 &&  error_residual_norm.u < param.newtonraphson.res_tol/* error_residual_norm.norm < prev_error_norm_m*/){  //TO DO
            std::cout<<"Converged"<<std::endl;
            print_footer();
	    break;
        }
//        assemble_system(param,newton_update);
//        make_constraints(new_iter,current_time_m);

        const std::pair<unsigned int,double>
                lin_solver_output = solve_linear_sys(param,newton_update);

        solution_delta += newton_update;

        std::cout << " | " << std::fixed << std::setprecision(3) << std::setw(7)
                            << std::scientific << lin_solver_output.first << "  "
                            << lin_solver_output.second << "  " << error_residual_norm.norm
                            << "  " << error_residual_norm.u << "  " << error_residual_norm.d
			    << "  " << std::endl;

    }
    AssertThrow (new_iter < param.newtonraphson.max_new_ite,
                   ExcMessage("No convergence in nonlinear solver!"));

}

template <int dim>
std::pair<unsigned int,double> Phasefield<dim>::solve_linear_sys(const AllParameters &param,BlockVector<double> &newton_update)
{
/*
  std::cout<<std::endl;
  std::cout<<"Initial Newton_update.block(d_dof_m): "<<newton_update.block(d_dof_m).l2_norm()<<std::endl;
  std::cout<<"Initial Newton_update.block(u_dof_m): "<<newton_update.block(u_dof_m).l2_norm()<<std::endl;
*/

  unsigned int lin_ite = 0;
  double lin_res = 0;
//  newton_update = 0;    // (solution)
  Vector<double> tmp1(newton_update.block(u_dof_m).size());
  Vector<double> tmp2(newton_update.block(u_dof_m).size());


  SolverControl           solver_control (tangent_matrix_m.block(d_dof_m,d_dof_m).m(),
                                          param.linearsolver.cg_tol*system_rhs_m.block(d_dof_m).l2_norm());
  GrowingVectorMemory<Vector<double> > GVM;
  SolverCG<Vector<double>>    cg (solver_control,GVM);
  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(tangent_matrix_m.block(d_dof_m,d_dof_m), param.linearsolver.relax_prm);
  cg.solve (tangent_matrix_m.block(d_dof_m,d_dof_m), newton_update.block(d_dof_m), system_rhs_m.block(d_dof_m),preconditioner);

/*
  SparseDirectUMFPACK k_dd;
  k_dd.initialize(tangent_matrix_m.block(d_dof_m,d_dof_m)); 
  k_dd.vmult(newton_update.block(d_dof_m),system_rhs_m.block(d_dof_m));
*/
  lin_ite = solver_control.last_step();
  lin_res = solver_control.last_value();

  constraints_m.distribute (newton_update);
/*
  std::cout<<std::endl;
  std::cout<<"Solved Newton_update.block(d_dof_m): "<<newton_update.block(d_dof_m).l2_norm()<<std::endl;
  std::cout<<"Initial Newton_update.block(u_dof_m): "<<newton_update.block(u_dof_m).l2_norm()<<std::endl;
*/
  SparseDirectUMFPACK k_uu;
  k_uu.initialize(tangent_matrix_m.block(u_dof_m,u_dof_m)); 

  tmp1 = system_rhs_m.block(u_dof_m);
//  std::cout<<"tmp1 [f_u]: "<<tmp1.l2_norm()<<std::endl;

  tangent_matrix_m.block(u_dof_m,d_dof_m).vmult(tmp2,newton_update.block(d_dof_m));
//  std::cout<<"tmp2 [k_ud*delta_d]: "<<tmp2.l2_norm()<<std::endl;

  tmp2 *= -1;
  tmp1 += tmp2;
//  std::cout<<"tmp1 [f_u - k_ud*delta_d]: "<<tmp1.l2_norm()<<std::endl;

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
void Phasefield<dim>::refine_grid (const AllParameters &param)
{
  Vector<float> estimated_error_per_cell (triangulation_m.n_active_cells());
  KellyErrorEstimator<dim>::estimate (dof_handler_m,
                                      QGauss<dim-1>(param.fesys.quad_order),
                                      typename FunctionMap<dim>::type(),
                                      solution_m,
                                      estimated_error_per_cell);
  GridRefinement::refine_and_coarsen_fixed_number (triangulation_m,
                                                   estimated_error_per_cell,
                                                   param.geometrymodel.act_ref, param.geometrymodel.act_cors);
 
  SolutionTransfer<dim,BlockVector<double>> solution_trans(dof_handler_m);
  BlockVector<double> tmp;
  tmp = old_solution_m;
  triangulation_m.prepare_coarsening_and_refinement();
  solution_trans.prepare_for_coarsening_and_refinement(tmp);

  triangulation_m.execute_coarsening_and_refinement ();

  setup_system();

  solution_trans.interpolate(tmp,old_solution_m);
  constraints_m.distribute(old_solution_m);

}


template <int dim>
void Phasefield<dim>::output_results (const double cycle) const
{

  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler_m);
  std::vector<std::string> solution_names;
  switch (dim)
    {
    case 1:
      solution_names.emplace_back("displacement");
      solution_names.emplace_back("phase_field");
      break;
    case 2:
      solution_names.emplace_back("x_displacement");
      solution_names.emplace_back("y_displacement");
      solution_names.emplace_back("phase_field");
      break;
    case 3:
      solution_names.emplace_back("x_displacement");
      solution_names.emplace_back("y_displacement");
      solution_names.emplace_back("z_displacement");
      solution_names.emplace_back("phase_field");
      break;
    default:
      Assert (false, ExcNotImplemented());
    }
  data_out.add_data_vector (solution_m, solution_names);
  data_out.build_patches ();
  std::ofstream output ("solution-" + std::to_string(cycle) + ".vtk");
  data_out.write_vtk (output);
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
    const unsigned int l_width = 80;
        for (unsigned int i = 0; i < l_width; ++i)
        {
            std::cout << "_";
        }
        std::cout << std::endl;

        std::cout << " SOLVER STEP"
                  << "|  LIN_IT   LIN_RES    RES_NORM    "
		  << " RES_U    RES_D    "                    
		  << std::endl;

        for (unsigned int i = 0; i < l_width; ++i)
        {
            std::cout << "_";
        }
        std::cout << std::endl;

}

template <int dim>
void Phasefield<dim>::print_footer(){
    const unsigned int l_width = 80;
        for (unsigned int i = 0; i < l_width; ++i)
        {
            std::cout << "_";
        }
        std::cout << std::endl;

        std::cout << "Errors:" << std::endl
                  << "Rhs_u: \t\t" << error_residual_norm.u << std::endl
                  << "Rhs_d: \t\t" << error_residual_norm.d << std::endl
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
void Phasefield<dim>::make_constraints(unsigned int &itr,const double time){
    std::cout<<" CST "<<std::flush;

    if(itr>=1){
        return;
    }

   constraints_m.clear();
   DoFTools::make_hanging_node_constraints (dof_handler_m,
                                           constraints_m);
/*
   VectorTools::interpolate_boundary_values (dof_handler_m,
                                              0,
                                              Functions::ZeroFunction<dim>(n_components_m),
                                              constraints_m
                                              ,fe_m.component_mask(u_extractor));
*/	
   std::vector<bool> component_mask(dim+1, false);
      
   //Tension test
   component_mask[0] = false;
   component_mask[1] = true;
   VectorTools::interpolate_boundary_values(dof_handler_m, 2,
                                               ZeroFunction<dim>(dim+1), constraints_m, component_mask);

   component_mask[0] = true;
   component_mask[1] = true;
   VectorTools::interpolate_boundary_values(dof_handler_m, 3,
                                               BoundaryTension<dim>(time), constraints_m, component_mask);
   constraints_m.close ();

   constraints_m.distribute(old_solution_m);   

}


template <int dim>
void Phasefield<dim>::run(const AllParameters &param){


    using namespace constants;

    BlockVector<double>       solution_delta(dofs_per_block_m);

    long double present_time_tol;
    present_time_tol = time_tol_fac * param.time.delta_t;
    current_time_m = param.time.start_time;
    unsigned int current_timestep = 0;

    std::cout<<"FESystem:n_blocks:"<<fe_m.n_blocks()<<std::endl;
    std::cout<<"FESystem:n_components:"<<fe_m.n_components()<<std::endl;
    std::cout<<"FESystem:n_base_elements:"<<fe_m.n_base_elements()<<std::endl;

    output_results(/*current_time_m*/current_timestep);
    current_time_m += param.time.delta_t;

    while (/*current_time_m < param.time.end_time + present_time_tol*/
	current_timestep < param.time.max_timesteps)
    {
        if(current_time_m >param.time.delta_t)
	{
            	refine_grid(param);
		solution_delta.reinit(dofs_per_block_m);
	}
/////////Printing////////
	std::cout<<"Parameter lambda: "<<param.materialmodel.lambda<<std::endl;
	std::cout<<"Parameter mu: "<<param.materialmodel.mu<<std::endl;
	std::cout<<"Parameter vis: "<<param.materialmodel.viscosity<<std::endl;
	std::cout<<"Parameter g_c: "<<param.pf.g_c<<std::endl;
	std::cout<<"Parameter l: "<<param.pf.l<<std::endl;
	std::cout<<"Parameter k: "<<param.pf.k<<std::endl;
	std::cout<<"Parameter delta_t: "<<param.time.delta_t<<std::endl;
///////////////////////
        
	solution_delta = 0.0;
        solve_nonlinear_newton(param,solution_delta);
        solution_m += solution_delta;				
        output_results(/*current_time_m*/current_timestep);
	
	old_old_history_m = old_history_m;
	old_history_m = std::max(old_old_history_m,/*get_history(param.materialmodel.lambda,param.materialmodel.mu,solution_m)*/-1.0);
	
        std::cout<<" solution.norm():"<<solution_m.l2_norm()<<std::endl;
	std::cout<<std::endl;
/*
        if(current_time_m >param.time.delta_t)
	{
            	refine_grid(param);
		setup_system();
		solution_delta.reinit(dofs_per_block_m);
	}
*/

	old_solution_m = solution_m;

        current_time_m += param.time.delta_t;
	++current_timestep;
//	compute_load(param.materialmodel.lambda,param.materialmodel.mu,solution_m);
   }

}

template class thesis::Phasefield<2>;
template class thesis::Phasefield<3>;

