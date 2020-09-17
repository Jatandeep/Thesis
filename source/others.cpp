#include "../include/others.h"
#include "../include/utilities.h"

using namespace dealii;
using namespace thesis;
using namespace parameters;

template <int dim>
double BoundaryTension<dim>::value (const Point<dim>  &,
                                    const unsigned int ) const
{}

template <int dim>
void BoundaryTension<dim>::vector_value (const Point<dim> &,
                                           Vector<double>   &vectorValue) const
{
  if(itr_ >0){
  vectorValue = 0.0;
//   vectorValue[1] = 0.0;
  }
  else
  {
	vectorValue[0] = 0.0;
	vectorValue[1] = load_ratio_ *u_total_;

  }
  
} 

template <int dim>
double BoundaryShear<dim>::value (const Point<dim>  &,
                                    const unsigned int) const
{}

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

template <int dim>
double InitialCrack<dim>::value (const Point<dim>  &p,
                                    const unsigned int component) const
{
// 0 = no crack
// 1 = crack 

if (component == dim)
   {
     if(itr_ >0)          //required
	    return 0.0;
    
    return 1;
   }    
}

template <int dim>
void InitialCrack<dim>::vector_value (const Point<dim> &p,
                                            Vector<double>   &values) const
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values (c) = InitialCrack<dim>::value (p, c);
}

template <int dim>
double InitialValuesCrack<dim>::value (const Point<dim>  &p,
                                    const unsigned int component) const
{
// 0 = no crack
// 1 = crack 
	//std::cout<<"min cell dia:"<<_min_cell_dia<<std::endl;
if (component == dim)
   {
     if ( (p(0) >= 0) && (p(0) <= 0.5 /*+ _min_cell_dia*/) && (p(1) >= 0.5-.01/*_min_cell_dia*/) && (p(1) <= 0.5+.01/*_min_cell_dia*/) )
        return 1.0;
      else
        return 0.0;
   }
return 0;
  
}

template <int dim>
void InitialValuesCrack<dim>::vector_value (const Point<dim> &p,
                                            Vector<double>   &values) const
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values (c) = InitialValuesCrack<dim>::value (p, c);
}


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
        
        if((support_points[idx][0] <= 0.5) 
        && (support_points[idx][0] >= 0.0) 
        && (support_points[idx][1] <= 0.5 + (min_cell_dia/10)) 
        && (support_points[idx][1] >= 0.5-(min_cell_dia/10)) )
        {
          global_index_m.push_back(idx);
        }
      }
    }
  }

/*
  for (const auto &cell : dof_handler_m.active_cell_iterators())
  {
    for (unsigned int vertex = 0;vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex)
    {
                Tensor<1, dim> cell_vertex = (cell->vertex(vertex));
                if (cell_vertex[0] <= 0.5 && cell_vertex[0] >= 0.0 
                    && cell_vertex[1] <= 0.5 + (min_cell_dia/10) && cell_vertex[1] >= 0.5-(min_cell_dia/10))
                {
                const unsigned int idx = cell->vertex_dof_index(vertex,2);
                global_index_m.push_back(idx);
                }
    }
  }
 */
}

template<int dim>
void Phasefield<dim>::get_constrained_initial_d(unsigned int itr_,const AllParameters &param)
{
  std::vector<unsigned int> local_dof_indices(fe_m.dofs_per_cell);

  constraints_m.clear();

  DoFTools::make_hanging_node_constraints (dof_handler_m,
                                            constraints_m);


  for(auto idx:global_index_m)
  {
   constraints_m.add_line(idx);
   if(current_time_m==param.time.delta_t)
    {
     if(itr_ >0)          
      {
       constraints_m.set_inhomogeneity(idx,0);
      }
     else
      {
       constraints_m.set_inhomogeneity(idx,1);
      }
    }
    else
    {
     constraints_m.set_inhomogeneity(idx,0);
    }
  }   
}

template <int dim>
double Reference_solution<dim>::value (const Point<dim>  &,
                                    const unsigned int ) const
{}

template <int dim>
void Reference_solution<dim>::vector_value (const Point<dim> &p,
                                           Vector<double>   &vectorValue) const
{
  using numbers::PI;
 
  const double x = p(0);
  const double y = p(1);

  const double k = 3-4*get_youngsM_poissonR(lambda_,mu_).second;
  const double K_I = load_ratio_*u_total_*sqrt(PI*0.5); 

  const double r = std::sqrt(x*x +y*y);
  const double theta = atan(y/x);
  
  if(itr_ >0){
  vectorValue = 0.0;
  }
  else
  {
	vectorValue[0] = (K_I/(4*mu_))*std::sqrt(r/(2*PI))*(-(2*k - 1)*sin(theta/2) - sin(3*theta/2));
	vectorValue[1] = (K_I/(4*mu_))*std::sqrt(r/(2*PI))*( (2*k + 1)*cos(theta/2) + cos(3*theta/2));
  }
  
}

template <int dim>
double Reference_solution_neg<dim>::value (const Point<dim>  &,
                                    const unsigned int ) const
{}

template <int dim>
void Reference_solution_neg<dim>::vector_value (const Point<dim> &p,
                                           Vector<double>   &vectorValue) const
{
  using numbers::PI;
 
  const double x = p(0);
  const double y = p(1);

  const double k = 3-4*get_youngsM_poissonR(lambda_,mu_).second;
  const double K_I = load_ratio_*u_total_*sqrt(PI*0.5); 

  const double r = std::sqrt(x*x +y*y);
  const double theta = atan(y/x);
  
  if(itr_ >0){
  vectorValue = 0.0;
  }
  else
  {
	vectorValue[0] = (K_I/(4*mu_))*std::sqrt(r/(2*PI))*(-(2*k - 1)*sin(theta/2) - sin(3*theta/2));
	vectorValue[1] = -(K_I/(4*mu_))*std::sqrt(r/(2*PI))*( (2*k + 1)*cos(theta/2) + cos(3*theta/2));
  }
  
}
/*
template <int dim>
double Reference_solution<dim>::value(const Point<dim> &   p,
                          const unsigned int component) const
{
    Assert(component <= 2 + 1, ExcIndexRange(component, 0, 2 + 1));
 
    using numbers::PI;
 
    const double x = p(0);
    const double y = p(1);
    
    const double K_I = 1e+4;
    const double mu = 87500;
    const double k = 2.20;

    //make get_theta function to ensure it is between 0 and 2 PI (like 215 degrees) and not between 0 and PI/2(or -PI/2 to PI/2)
    
 
    const double r = std::sqrt(x*x +y*y);
    const double theta = atan(y/x);

    if (component == 0)
      return ( (K_I/(4*mu))*std::sqrt(r/(2*PI))*(-(2*k - 1)*sin(theta/2) - sin(3*theta/2)) ); 
    if (component == 1)
      return ( (K_I/(4*mu))*std::sqrt(r/(2*PI))*( (2*k + 1)*cos(theta/2) + cos(3*theta/2)) );
   
    return 0;
}
*/
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

double get_deg_func(const double d)
{
return ( pow((1-d),2.0) );
}

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
  
  for (const auto &cell : dof_handler_m.active_cell_iterators())
    {

        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
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
                  /*
                  if (cell->face(face)->boundary_id() == 0)
                  {
                    Point<dim> point_1;
                    point_1 = cell->face(face)->barycenter();
                    std::cout<<"Point:barycenter:  "<<point_1[0]<<"  "<<point_1[1]<<std::endl;

                  }*/
	        }
    }
    
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

template <int dim>
void ElasticBodyForce<dim>::vector_value (const Point<dim> &points,
                                  Vector<double>  &values) const
{
    Point<dim> point_1, point_2;
    point_1(0) = 0.5;
    point_2(0) = -0.5;

    if (((points-point_1).norm_square() < 0.2*0.2) ||
            ((points-point_2).norm_square() < 0.2*0.2))
       values[0] = 1.0;
    else
       values[0] = 0.0;
    if (points.norm_square() < 0.2*0.2)
       values[1] = 1.0;
    else
       values[1] = 0.0;

}

double get_stress_intensity_const(const std::string test)
{
  using numbers::PI;

if(test == "tension")
   {
     double f_a_b = 1.12 - 0.23*(0.5) + 10.56*pow(0.5,2.0) - 21.74*pow(0.5,3.0) + 30.42*pow(0.5,4.0) ;
           
     std::cout<<"F(a/b):"<<f_a_b<<std::endl;
     return ( std::sqrt(PI*0.5)*f_a_b );
            
   }
else if(test == "shear")
  {
    return 0;
  }
}


std::pair<double,double> get_youngsM_poissonR(const double lambda,const double mu)
{
  double YoungsM,PoissonR;
 
  PoissonR = ((3*lambda - 2*mu)/(6*lambda + 2*mu) );
  YoungsM = 2*mu*(1+PoissonR);

  return std::make_pair(YoungsM,PoissonR);
}
/*
template<int dim>
double get_critical_load(Tensor<2,dim> stress)
{
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
                  if (cell->face(face)->at_boundary() && cell->face(face)->boundary_id() == 4)
                  {
                        fe_values_face.reinit(cell, face);
        					                                                                                  
                        for (unsigned int q_point = 0; q_point < n_face_q_points;++q_point)
                        {
                          load_value += stress*fe_values_face.normal_vector(q_point)*fe_values_face.JxW(q_point);
                        }
                  }
	        }
    }
  return load_value;
}
*/


template class thesis::ElasticBodyForce<2>;
template class thesis::ElasticBodyForce<3>;
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
template class thesis::InitialValuesCrack<2>;
template class thesis::InitialValuesCrack<3>;
template class thesis::Reference_solution<2>;
template class thesis::Reference_solution<3>;
template class thesis::Reference_solution_neg<2>;
template class thesis::Reference_solution_neg<3>;
template void thesis::Phasefield<2>::compute_load(const AllParameters &,BlockVector<double> &);
template void thesis::Phasefield<3>::compute_load(const AllParameters &,BlockVector<double> &);
template void thesis::Phasefield<2>::get_constrained_initial_d(unsigned int,const AllParameters &);
template void thesis::Phasefield<3>::get_constrained_initial_d(unsigned int,const AllParameters &);
template void thesis::Phasefield<2>::extract_initialcrack_d_index(const double,const AllParameters &);                                                              
template void thesis::Phasefield<3>::extract_initialcrack_d_index(const double,const AllParameters &);
                                                              