#pragma once
#include <iostream>
#include <deal.II/lac/vector.h>
#include <deal.II/base/point.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/fe/fe_values.h>

namespace thesis
{
    template <int dim>
    void right_hand_side (const std::vector<dealii::Point<dim> > &points,
                          std::vector<dealii::Tensor<1, dim> >   &values);
    template <int dim>
    inline dealii::SymmetricTensor<2, dim> get_strain(const dealii::FEValues<dim> &fe_values,
                                              const unsigned int   shape_func,
                                              const unsigned int   q_point);

};









