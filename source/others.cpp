#include "../include/ElasticProblem.h"

using namespace dealii;
using namespace thesis;

template <int dim>
SymmetricTensor<4, dim> get_stress_strain_tensor(const double lambda,
                                                 const double mu)
{
  SymmetricTensor<4, dim> tmp;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        for (unsigned int l = 0; l < dim; ++l)
          tmp[i][j][k][l] = (((i == j) && (k == l)) ? lambda:0)
                            + (((i == k) && (j == l)) ? mu:0)
                            + (((i == l) && (j == k)) ? mu:0);
  return tmp;
}

template <int dim>
inline SymmetricTensor<2, dim> get_strain(const FEValues<dim> &fe_values,
                                          const unsigned int   shape_func,
                                          const unsigned int   q_point)
{
  SymmetricTensor<2, dim> temp;

  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i ; j < dim; ++j)
      temp[i][j] =
        (fe_values.shape_grad_component(shape_func, q_point, i)[j] +
         fe_values.shape_grad_component(shape_func, q_point, j)[i]) /
        2;
  return temp;
}

template <int dim>
const SymmetricTensor<4, dim> ElasticProblem<dim>::stress_strain_tensor =
        get_stress_strain_tensor<dim>(1, 1);


template <int dim>
void right_hand_side (const std::vector<Point<dim> > &points,
                      std::vector<Tensor<1, dim> >   &values)
{
  Assert (values.size() == points.size(),
          ExcDimensionMismatch (values.size(), points.size()));
  Assert (dim >= 2, ExcNotImplemented());
  Point<dim> point_1, point_2;
  point_1(0) = 0.5;
  point_2(0) = -0.5;
  for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
    {
      if (((points[point_n]-point_1).norm_square() < 0.2*0.2) ||
          ((points[point_n]-point_2).norm_square() < 0.2*0.2))
        values[point_n][0] = 1.0;
      else
        values[point_n][0] = 0.0;
      if (points[point_n].norm_square() < 0.2*0.2)
        values[point_n][1] = 1.0;
      else
        values[point_n][1] = 0.0;
    }
}





template class thesis::ElasticProblem<2>;
