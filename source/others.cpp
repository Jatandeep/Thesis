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
/*
template <int dim>
struct InnerPreconditioner;

template <>
class InnerPreconditioner<2>
{
  using type = SparseDirectUMFPACK;
};

template <>
class InnerPreconditioner<3>
{
  using type = SparseILU<double>;
};
*/

template <class MatrixType, class PreconditionerType>
InverseMatrix<MatrixType, PreconditionerType>::InverseMatrix(const MatrixType &m,
    							     const PreconditionerType &preconditioner)
    : matrix(&m)
    , preconditioner(&preconditioner)
  {}

template <class MatrixType, class PreconditionerType>
void InverseMatrix<MatrixType, PreconditionerType>::vmult (Vector<double> &dst,
    							   const Vector<double> &src) const
  {
    SolverControl solver_control(src.size(), 1e-6 * src.l2_norm());
    SolverCG<>    cg(solver_control);
    dst = 0;
    cg.solve(*matrix, dst, src, *preconditioner);
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
double get_history(const double lambda
	  ,const double mu
	  ,SymmetricTensor<2,dim> &eps)
{
	
	double history;
	double tr_eps = trace(eps);
	double tr_eps_plus = (tr_eps>0) ? tr_eps:0;
	double tr_epsplus_sq = get_epsplus_sq(eps); 	


	history = 0.5*lambda*tr_eps_plus*tr_eps_plus + mu*tr_epsplus_sq;

return history;
}




template class thesis::BoundaryForce<2>;
template class thesis::BoundaryForce<3>;
template class thesis::BoundaryTension<2>;
template class thesis::BoundaryTension<3>;
//template struct thesis::InnerPreconditioner<2>;
//template struct thesis::InnerPreconditioner<3>;
template double get_history(const double,const double,SymmetricTensor<2,2>&);
template double get_history(const double,const double,SymmetricTensor<2,3>&);

