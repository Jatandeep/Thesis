#include "../include/ElasticProblem.h"
#include "../include/others.h"
#include "../include/Phasefield.h"
using namespace dealii;
using namespace thesis;

template <int dim>
double BoundaryTension<dim>::value (const Point<dim>  &p,
                                    const unsigned int component) const
{
  Assert (component < this->n_components,
          ExcIndexRange (component, 0, this->n_components));

  double u_step_per_timestep = 1.0;

  if (component == 1)
    {
      return ( ((p(1) == 1.0) && (p(0) <= 1.0) && (p(0) >= 0.0))
               ? (time_ * u_step_per_timestep) : 0 );

    }

  return 0;
}

template <int dim>
void BoundaryTension<dim>::vector_value (const Point<dim> &p,
                                           Vector<double>   &values) const
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values (c) = BoundaryTension<dim>::value (p, c);
}


template <int dim>
BoundaryForce<dim>::BoundaryForce ()
:Function<dim>(dim+1){}
   
template <int dim>
double BoundaryForce<dim>::value (const Point<dim> &,
                               const unsigned int) const
{
   return 0.;
}

template <int dim>
void BoundaryForce<dim>::vector_value (const Point<dim> &p,
                                      Vector<double> &values) const
{
      for (unsigned int c = 0; c < this->n_components; ++c)
        values(c) = BoundaryForce<dim>::value(p, c);
}


template <int dim>
double get_epsplus_sq(SymmetricTensor<2,dim> &eps)
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
double Phasefield<dim>::get_history(const double lambda
	  	  ,const double mu
	  	  ,BlockVector<double> &solution)
{
	
	FEValues<dim> fe_values (fe_m, quadrature_formula_m,
                           update_values   | update_gradients |
                           update_quadrature_points | update_JxW_values);
	const unsigned int   n_q_points    = quadrature_formula_m.size();
 	double history;

	for (const auto &cell : dof_handler_m.active_cell_iterators())
	{
		fe_values.reinit(cell);
		
		std::vector<SymmetricTensor<2,dim>> eps(n_q_points);
		fe_values[u_extractor].get_function_symmetric_gradients(solution,eps);

        	for (unsigned int q = 0; q < n_q_points; ++q)
		{
			double tr_eps = trace(eps[q]);
			double tr_eps_plus = (tr_eps>0) ? tr_eps:0;
			double tr_epsplus_sq = get_epsplus_sq(eps[q]); 	

			history += 0.5*lambda*tr_eps_plus*tr_eps_plus + mu*tr_epsplus_sq;	//+= ??
		}
	}
return history;
}


template <int dim>
void Phasefield<dim>::compute_load(const double lambda
                 ,const double mu
                 ,dealii::BlockVector<double> &solution)
{

  FEFaceValues<dim> fe_values_face(fe_m, face_quadrature_formula_m,
                                   update_values | update_quadrature_points |
                                   update_normal_vectors | update_JxW_values);

  const unsigned int   dofs_per_cell = fe_m.dofs_per_cell;
  const unsigned int   n_face_q_points    = face_quadrature_formula_m.size();

  Tensor<1,dim> load_value;
  Tensor<1,dim> load_value_1;

  for (const auto &cell : dof_handler_m.active_cell_iterators())
    {

        fe_values_face.reinit(cell);

        std::vector<SymmetricTensor<2,dim>> epsilon_vals(n_face_q_points);
        fe_values_face[u_extractor].get_function_symmetric_gradients(solution,epsilon_vals);

	for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
        {
                  if (cell->face(face)->at_boundary() && cell->face(face)->boundary_id() == 3)
                  {
                        fe_values_face.reinit(cell, face);
        
                        for (unsigned int q_point = 0; q_point < n_face_q_points;++q_point)
                        {
              
                                Tensor<2,dim> stress_display;
                                double tr_eps;
                                tr_eps = trace(epsilon_vals[q_point]);
                                stress_display = lambda*tr_eps*unit_symmetric_tensor<dim>()
                                                + 2*mu*epsilon_vals[q_point];
                                load_value += stress_display*fe_values_face.normal_vector(q_point)*fe_values_face.JxW(q_point);

                        }
                  }
		 if(cell->face(face)->at_boundary() && (/*cell->face(face)->boundary_id() == 0 ||
                                                        cell->face(face)->boundary_id() == 1 ||*/
                                                        cell->face(face)->boundary_id() == 2) )
                 {
                        fe_values_face.reinit(cell, face);
                                for (unsigned int q_point = 0; q_point < n_face_q_points;++q_point)
                        {

                                Tensor<2,dim> stress_display_1;
                                double tr_eps;
                                tr_eps = trace(epsilon_vals[q_point]);
                                stress_display_1 = lambda*tr_eps*unit_symmetric_tensor<dim>()
                                                + 2*mu*epsilon_vals[q_point];
                                load_value_1 += stress_display_1*fe_values_face.normal_vector(q_point)*fe_values_face.JxW(q_point);
                        }
                }
        }
  }
/*
std::cout<<std::endl;
std::cout<<"Boundary3_Load_x: "<<load_value[0]<<std::endl;
std::cout<<"Boundary3_Load_y: "<<load_value[1]<<std::endl;
std::cout<<"Boundary012_Load_x: "<<load_value_1[0]<<std::endl;
std::cout<<"Boundary012_Load_y: "<<load_value_1[1]<<std::endl;
*/
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

template class thesis::ElasticBodyForce<2>;
template class thesis::ElasticBodyForce<3>;
template class thesis::BoundaryForce<2>;
template class thesis::BoundaryForce<3>;
template class thesis::BoundaryTension<2>;
template class thesis::BoundaryTension<3>;
template double thesis::Phasefield<2>::get_history(const double,const double,BlockVector<double>&);
template double thesis::Phasefield<3>::get_history(const double,const double,BlockVector<double>&);
//template double get_history(const double,const double,SymmetricTensor<2,3>&);
//template void compute_load(const double,const double,BlockVector<double>&);

