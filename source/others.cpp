#include "../include/others.h"
#include "../include/utilities.h"

using namespace dealii;
using namespace thesis;
using namespace parameters;

/*A dealii type function for Tension test for prescribing load on top boundary*/
template <int dim>
void BoundaryTension<dim>::vector_value (const Point<dim> &,
                                           Vector<double>   &vectorValue) const
{
  if(itr_ >0){
  vectorValue = 0.0;
  }
  else
  {
	vectorValue[0] = 0.0;
	vectorValue[1] = load_ratio_ *u_total_;
  }  
} 

/*A dealii type function for shear test for prescribing load on top boundary*/
template <int dim>
void BoundaryShear<dim>::vector_value (const Point<dim> &,
                                           Vector<double>   &vectorValue) const
{
  if(itr_ >0)
	  vectorValue = 0.0;
  else
  {
	vectorValue[0] = (1)*load_ratio_ *u_total_;
	vectorValue[1] = 0.0;
  }
  
} 

/*For M_Id, this dealii type function assigns d=1 on crack boundary*/
template <int dim>
double InitialCrack<dim>::value (const Point<dim>  &,
                                    const unsigned int component) const
{
/* 0 = no crack; 1 = crack*/ 

if (component == dim)
   {
     if(itr_ >0)          
	    return 0.0;
    
    return 1;
   }  
  return 0;   
}

template <int dim>
void InitialCrack<dim>::vector_value (const Point<dim> &p,
                                            Vector<double>   &values) const
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values (c) = InitialCrack<dim>::value (p, c);
}

/*For P_I, this function extracts global node ids from cells where crack is to be placed*/   
template<int dim>
void Phasefield<dim>::extract_initialcrack_d_index(const double min_cell_dia,const AllParameters &param)                                                              
{

  constraints_m.clear();
  DoFTools::make_hanging_node_constraints (dof_handler_m,
                                            constraints_m);
                                            
  std::vector<Point<dim>> support_points(dof_handler_m.n_dofs());
 
  MappingQ<dim> mapping(param.fesys.fe_degree,true);
 
  DoFTools::map_dofs_to_support_points<dim>(mapping
                                           ,dof_handler_m
                                           ,support_points);
  
  std::vector<unsigned int> local_dof_indices(fe_m.dofs_per_cell);
  for (const auto &cell : dof_handler_m.active_cell_iterators())
  {
    cell->get_dof_indices(local_dof_indices);
    for (unsigned int i=0; i<fe_m.dofs_per_cell; ++i)
    {
      const unsigned int comp_i = fe_m.system_to_component_index(i).first;

      if (comp_i != dim)
      continue; // only look at phase field
      else
      {
        const unsigned int idx = local_dof_indices[i];
        
        if(param.mod_strategy.comp_strategy=="benchmarks")
        {
          if(param.mod_strategy.pi_strategy=="elements")
          {
            if((support_points[idx][0] <= param.geometrymodel.b) 
            && (support_points[idx][0] >= 0.0) 
            && (support_points[idx][1] <= (param.geometrymodel.a/2) + min_cell_dia) 
            && (support_points[idx][1] >= (param.geometrymodel.a/2) - min_cell_dia) )
            {
              global_index_m.push_back(idx);
            }
          }
          if(param.mod_strategy.pi_strategy=="nodes")
          {
            if((support_points[idx][0] <= param.geometrymodel.b) 
            && (support_points[idx][0] >= 0.0) 
            && (support_points[idx][1] <= (param.geometrymodel.a/2) + min_cell_dia/10) 
            && (support_points[idx][1] >= (param.geometrymodel.a/2) - min_cell_dia/10) )
            {
              global_index_m.push_back(idx);
            }
          }
        }
       
        if(param.mod_strategy.comp_strategy=="lefm_mode_I")
        {
          if((support_points[idx][0] <= 0.0) 
          && (support_points[idx][0] >= -(param.geometrymodel.a/2)) 
          && (support_points[idx][1] <= min_cell_dia) 
          && (support_points[idx][1] >= -min_cell_dia) )
          {
            global_index_m.push_back(idx);
          }
        }
         
      }
    }
  }
  std::sort(global_index_m.begin(),global_index_m.end());
  global_index_m.erase(std::unique(global_index_m.begin(),global_index_m.end())
                      ,global_index_m.end());
  std::cout<<"SIze global vector"<<global_index_m.size()<<std::endl;

}

/*For P_I method, this function prescribe d=1 values on selected nodes for initial time step and maintains them throughout simulation*/
template<int dim>
void Phasefield<dim>::get_constrained_initial_d(unsigned int itr_,const AllParameters &param)
{
  std::vector<unsigned int> local_dof_indices(fe_m.dofs_per_cell);

  constraints_m.clear();

  DoFTools::make_hanging_node_constraints (dof_handler_m,
                                            constraints_m);
  
   
  if(current_time_m==param.time.delta_t)
  {
      for(auto idx:global_index_m)
      {
        constraints_m.add_line(idx);
      if(itr_ >0)          
        {
        constraints_m.set_inhomogeneity(idx,0);
        }
      else
        {
        constraints_m.set_inhomogeneity(idx,1);
        }
      }
      
  }
  else
  {
    for(auto idx:global_index_m)
    {
     constraints_m.add_line(idx);
     constraints_m.set_inhomogeneity(idx,0);
    }

  }
     
}

/*For LEFM, this dealii type function prescribe the loading conditions on boundaries*/
template <int dim>
void Reference_solution<dim>::vector_value (const Point<dim> &p,
                                           Vector<double>   &vectorValue) const
{

  if(itr_ >0){
  vectorValue = 0.0;
  }
  else
  {
    using numbers::PI;
  
    const double x = p(0);
    const double y = p(1);      

    //Griffith crack growth criteria:
    const double k = 3-4*get_youngsM_poissonR(lambda_,mu_).second;
    const double fracture_toughness = std::sqrt((g_c_*get_youngsM_poissonR(lambda_,mu_).first)
                                                /(1-std::pow(get_youngsM_poissonR(lambda_,mu_).second,2)));
    
    const double K_I = fracture_toughness/steps_;

    double r = std::sqrt(x*x +y*y);
    double theta_radian = atan2(y,x);

    vectorValue[0] = (K_I/(2*mu_))*std::sqrt(r/(2*PI))*(cos(theta_radian/2)*(k-cos(theta_radian)));
    vectorValue[1] = (K_I/(2*mu_))*std::sqrt(r/(2*PI))*(sin(theta_radian/2)*(k-cos(theta_radian)) );
  }
  
}

/*To be used in get_history() and get_energy_density_plus()*/
template <int dim>
double get_epsplus_sq(const SymmetricTensor<2,dim> &eps)
{
	std::array<double,std::integral_constant<int,dim>::value> eps_eigenval;
	eps_eigenval = eigenvalues(eps);

	std::array<double,std::integral_constant<int,dim>::value> eps_eigenval_plus;

	for(unsigned int i=0;i<dim;++i){
		eps_eigenval_plus[i]=(eps_eigenval[i]>0) ? eps_eigenval[i]:0;
	}

	double result=0;
	for(unsigned int i=0;i<dim;++i){
	result += eps_eigenval_plus[i]*eps_eigenval_plus[i];
	}
	
return result;
}

/*To be used in get_energy_density_minus()*/
template <int dim>
double get_epsminus_sq(SymmetricTensor<2,dim> &eps)
{
	std::array<double,std::integral_constant<int,dim>::value> eps_eigenval;
	eps_eigenval = eigenvalues(eps);

	std::array<double,std::integral_constant<int,dim>::value> eps_eigenval_minus;
	for(unsigned int i=0;i<dim;++i){
		eps_eigenval_minus[i]=(eps_eigenval[i]>0) ? 0:eps_eigenval[i];
	}

	double result=0;
	for(unsigned int i=0;i<dim;++i){
	result += eps_eigenval_minus[i]*eps_eigenval_minus[i];
	}
	
return result;
}

/*Data member function for calculating the history function for present time step*/
template <int dim>
double Phasefield<dim>::get_history(const double lambda
	  	  ,const double mu
		  ,const SymmetricTensor<2,dim> &eps)const
{
	double history = 0;
			
	double tr_eps = trace(eps);
	double tr_eps_plus = (tr_eps>0) ? tr_eps:0;
	double tr_epsplus_sq = get_epsplus_sq(eps); 	

	history = 0.5*lambda*tr_eps_plus*tr_eps_plus + mu*tr_epsplus_sq;	

return history;
}

/*Generate energy density plus values for calculating total Elastic energy(for statistics file)*/
template <int dim>
double get_energy_density_plus(const double lambda
			  	  	       ,const double mu
				  	       ,SymmetricTensor<2,dim> &eps)
{
	double energy_plus = 0;
			
	double tr_eps = trace(eps);
	double tr_eps_plus = (tr_eps>0) ? tr_eps:0;
	double tr_epsplus_sq = get_epsplus_sq(eps); 	

	energy_plus = 0.5*lambda*tr_eps_plus*tr_eps_plus + mu*tr_epsplus_sq;	


return energy_plus;
}

/*Generate energy density minus values for calculating total Elastic energy(for statistics file)*/
template <int dim>
double get_energy_density_minus(const double lambda
			  	  	       ,const double mu
				  	       ,SymmetricTensor<2,dim> &eps)
{
	double energy_minus = 0;
			
	double tr_eps = trace(eps);
	double tr_eps_minus = (tr_eps>0) ? 0:tr_eps;
	double tr_epsminus_sq = get_epsminus_sq(eps); 	

	energy_minus = 0.5*lambda*tr_eps_minus*tr_eps_minus + mu*tr_epsminus_sq;	

return energy_minus;
}

/*Define the degradation function to be used in formulation*/
double get_deg_func(const double d)
{
return ( pow((1-d),2.0) );
}

/*Generating load on boundaries for being written to statistics file*/
template <int dim>
void Phasefield<dim>::compute_load(const AllParameters &param
                 ,dealii::BlockVector<double> &solution)
{
  double lambda = param.materialmodel.lambda;
  double mu = param.materialmodel.mu;

  FEFaceValues<dim> fe_values_face(fe_m, face_quadrature_formula_m,
                                   update_values | update_quadrature_points |
                                   update_normal_vectors | update_gradients |
                        				   update_JxW_values);

  const unsigned int   n_face_q_points    = face_quadrature_formula_m.size();

  Tensor<1,dim> load_value;
  Tensor<1,dim> load_value_lefm_1,load_value_lefm_2,load_value_lefm_3,load_value_lefm_4;
  
  for (const auto &cell : dof_handler_m.active_cell_iterators())
    {
        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
        {
                  if(param.mod_strategy.comp_strategy=="benchmarks")
                  {
                    if (cell->face(face)->at_boundary() && cell->face(face)->boundary_id() == 4)
                    {
                          fe_values_face.reinit(cell, face);
              
                          std::vector<SymmetricTensor<2,dim>> epsilon_vals(n_face_q_points);
                          fe_values_face[u_extractor].get_function_symmetric_gradients(solution,epsilon_vals);

                          std::vector<double> d_vals(n_face_q_points);
                          fe_values_face[d_extractor].get_function_values(solution,d_vals);
                                                                    
                          for (unsigned int q_point = 0; q_point < n_face_q_points;++q_point)
                          {
                
                                  Tensor<2,dim> stress_display;
                                  
                                  stress_display = (get_deg_func(d_vals[q_point]) + param.pf.k)
                                                  *get_stress_plus(lambda,mu,epsilon_vals[q_point])
                                                  + get_stress_minus(lambda,mu,epsilon_vals[q_point]);
                                  load_value += stress_display*fe_values_face.normal_vector(q_point)*fe_values_face.JxW(q_point);

                          }
                    
                    }
                  }
                  if(param.mod_strategy.comp_strategy=="lefm_mode_I")
                  {
                    if (cell->face(face)->at_boundary() && cell->face(face)->boundary_id() == 1)
                    {
                          fe_values_face.reinit(cell, face);
              
                          std::vector<SymmetricTensor<2,dim>> epsilon_vals(n_face_q_points);
                          fe_values_face[u_extractor].get_function_symmetric_gradients(solution,epsilon_vals);

                          std::vector<double> d_vals(n_face_q_points);
                          fe_values_face[d_extractor].get_function_values(solution,d_vals);
                                                                    
                          for (unsigned int q_point = 0; q_point < n_face_q_points;++q_point)
                          {
                
                                  Tensor<2,dim> stress_display;
                                  
                                  stress_display = (get_deg_func(d_vals[q_point]) + param.pf.k)
                                                  *get_stress_plus(lambda,mu,epsilon_vals[q_point])
                                                  + get_stress_minus(lambda,mu,epsilon_vals[q_point]);
                                  load_value_lefm_1 += stress_display*fe_values_face.normal_vector(q_point)*fe_values_face.JxW(q_point);

                          }
                    
                    }
                    if (cell->face(face)->at_boundary() && cell->face(face)->boundary_id() == 2)
                    {
                          fe_values_face.reinit(cell, face);
              
                          std::vector<SymmetricTensor<2,dim>> epsilon_vals(n_face_q_points);
                          fe_values_face[u_extractor].get_function_symmetric_gradients(solution,epsilon_vals);

                          std::vector<double> d_vals(n_face_q_points);
                          fe_values_face[d_extractor].get_function_values(solution,d_vals);
                                                                    
                          for (unsigned int q_point = 0; q_point < n_face_q_points;++q_point)
                          {
                
                                  Tensor<2,dim> stress_display;
                                  
                                  stress_display = (get_deg_func(d_vals[q_point]) + param.pf.k)
                                                  *get_stress_plus(lambda,mu,epsilon_vals[q_point])
                                                  + get_stress_minus(lambda,mu,epsilon_vals[q_point]);
                                  load_value_lefm_2 += stress_display*fe_values_face.normal_vector(q_point)*fe_values_face.JxW(q_point);

                          }
                    
                    }
                    if (cell->face(face)->at_boundary() && cell->face(face)->boundary_id() == 3)
                    {
                          fe_values_face.reinit(cell, face);
              
                          std::vector<SymmetricTensor<2,dim>> epsilon_vals(n_face_q_points);
                          fe_values_face[u_extractor].get_function_symmetric_gradients(solution,epsilon_vals);

                          std::vector<double> d_vals(n_face_q_points);
                          fe_values_face[d_extractor].get_function_values(solution,d_vals);
                                                                    
                          for (unsigned int q_point = 0; q_point < n_face_q_points;++q_point)
                          {
                
                                  Tensor<2,dim> stress_display;
                                  
                                  stress_display = (get_deg_func(d_vals[q_point]) + param.pf.k)
                                                  *get_stress_plus(lambda,mu,epsilon_vals[q_point])
                                                  + get_stress_minus(lambda,mu,epsilon_vals[q_point]);
                                  load_value_lefm_3 += stress_display*fe_values_face.normal_vector(q_point)*fe_values_face.JxW(q_point);

                          }
                    
                    }
                    if (cell->face(face)->at_boundary() && cell->face(face)->boundary_id() == 4)
                    {
                          fe_values_face.reinit(cell, face);
              
                          std::vector<SymmetricTensor<2,dim>> epsilon_vals(n_face_q_points);
                          fe_values_face[u_extractor].get_function_symmetric_gradients(solution,epsilon_vals);

                          std::vector<double> d_vals(n_face_q_points);
                          fe_values_face[d_extractor].get_function_values(solution,d_vals);
                                                                    
                          for (unsigned int q_point = 0; q_point < n_face_q_points;++q_point)
                          {
                
                                  Tensor<2,dim> stress_display;
                                  
                                  stress_display = (get_deg_func(d_vals[q_point]) + param.pf.k)
                                                  *get_stress_plus(lambda,mu,epsilon_vals[q_point])
                                                  + get_stress_minus(lambda,mu,epsilon_vals[q_point]);
                                  load_value_lefm_4 += stress_display*fe_values_face.normal_vector(q_point)*fe_values_face.JxW(q_point);

                          }
                    
                    }
                  }
	        }
    }
  if(param.mod_strategy.comp_strategy=="benchmarks")
  {  
    if(param.test_case.test == "tension"){
      double load_y = load_value[1];
      statistics.add_value("Load y",load_y);
    }
    else if(param.test_case.test == "shear"){
      load_value *= 1.0;
      double load_x = load_value[0];
      statistics.add_value("Load x",load_x);
    }
  }
  if(param.mod_strategy.comp_strategy=="lefm_mode_I")
  {
    if(param.test_case.test == "tension"){
      double load_y_1 = load_value_lefm_1[1];
      statistics.add_value("Load y Left",load_y_1);

      double load_y_2 = load_value_lefm_2[1];
      statistics.add_value("Load y Right",load_y_2);
      
      double load_y_3 = load_value_lefm_3[1];
      statistics.add_value("Load y Bottom",load_y_3);
      
      double load_y_4 = load_value_lefm_4[1];
      statistics.add_value("Load y Top",load_y_4);

    }
  }
}

/*Additional functional to generate a visualization file when load becomes 0 for tension case*/
template <int dim>
double Phasefield<dim>::compute_end_load(const AllParameters &param
                 ,dealii::BlockVector<double> &solution)
{
  double lambda = param.materialmodel.lambda;
  double mu = param.materialmodel.mu;

  FEFaceValues<dim> fe_values_face(fe_m, face_quadrature_formula_m,
                                   update_values | update_quadrature_points |
                                   update_normal_vectors | update_gradients |
                        				   update_JxW_values);

  const unsigned int   n_face_q_points    = face_quadrature_formula_m.size();

  Tensor<1,dim> load_value;
  
  for (const auto &cell : dof_handler_m.active_cell_iterators())
    {
        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
        {
                  if(param.mod_strategy.comp_strategy=="benchmarks")
                  {
                    if (cell->face(face)->at_boundary() && cell->face(face)->boundary_id() == 4)
                    {
                          fe_values_face.reinit(cell, face);
              
                          std::vector<SymmetricTensor<2,dim>> epsilon_vals(n_face_q_points);
                          fe_values_face[u_extractor].get_function_symmetric_gradients(solution,epsilon_vals);

                          std::vector<double> d_vals(n_face_q_points);
                          fe_values_face[d_extractor].get_function_values(solution,d_vals);
                                                                    
                          for (unsigned int q_point = 0; q_point < n_face_q_points;++q_point)
                          {
                
                                  Tensor<2,dim> stress_display;
                                  
                                  stress_display = (get_deg_func(d_vals[q_point]) + param.pf.k)
                                                  *get_stress_plus(lambda,mu,epsilon_vals[q_point])
                                                  + get_stress_minus(lambda,mu,epsilon_vals[q_point]);
                                  load_value += stress_display*fe_values_face.normal_vector(q_point)*fe_values_face.JxW(q_point);

                          }
                    
                    }
                  }
        }
    }
      
    if(param.test_case.test == "tension"){
      return load_value[1];
    } 
    else
    return 1;
}

/*Generate Youngs modulus and poisson ratio from lame parameters*/
std::pair<double,double> get_youngsM_poissonR(const double lambda,const double mu)
{
  double YoungsM,PoissonR;
 
  PoissonR = ((3*lambda - 2*mu)/(6*lambda + 2*mu) );
  YoungsM = 2*mu*(1+PoissonR);

  return std::make_pair(YoungsM,PoissonR);
}


template class thesis::BoundaryTension<2>;
template class thesis::BoundaryTension<3>;
template class thesis::BoundaryShear<2>;
template class thesis::BoundaryShear<3>;
template double thesis::Phasefield<2>::get_history(const double,const double,const SymmetricTensor<2,2>&)const;
template double thesis::Phasefield<3>::get_history(const double,const double,const SymmetricTensor<2,3>&)const;
template double get_energy_density_plus(const double,const double,SymmetricTensor<2,2>&);
template double get_energy_density_plus(const double,const double,SymmetricTensor<2,3>&);
template double get_energy_density_minus(const double,const double,SymmetricTensor<2,2>&);
template double get_energy_density_minus(const double,const double,SymmetricTensor<2,3>&);
template class thesis::InitialCrack<2>;
template class thesis::InitialCrack<3>;
template class thesis::Reference_solution<2>;
template class thesis::Reference_solution<3>;
template void thesis::Phasefield<2>::compute_load(const AllParameters &,BlockVector<double> &);
template void thesis::Phasefield<3>::compute_load(const AllParameters &,BlockVector<double> &);
template double thesis::Phasefield<2>::compute_end_load(const AllParameters &,BlockVector<double> &);
template double thesis::Phasefield<3>::compute_end_load(const AllParameters &,BlockVector<double> &);
template void thesis::Phasefield<2>::get_constrained_initial_d(unsigned int,const AllParameters &);
template void thesis::Phasefield<3>::get_constrained_initial_d(unsigned int,const AllParameters &);
template void thesis::Phasefield<2>::extract_initialcrack_d_index(const double,const AllParameters &);                                                              
template void thesis::Phasefield<3>::extract_initialcrack_d_index(const double,const AllParameters &);
