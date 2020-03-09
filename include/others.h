#pragma once
#include <iostream>
#include <deal.II/lac/vector.h>
#include <deal.II/base/point.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>
#include "../include/ElasticProblem.h"

namespace thesis
{

    template <int dim>
    class Others:public dealii::Function<dim>{
    public:
        Others():dealii::Function<dim>(1){}

        virtual double value (const dealii::Point<dim> &p,unsigned int component =0)const;

        virtual void vector_value (const dealii::Point<dim> &p,
                             dealii::Vector<double>  &value/*,double &step_frac*/)const;

    };


//    template <int dim>
//    class Display{
//    public:


//    };
/*
    template <int dim>
    class Error{
    public:
        Error();
        double u;
        void reset();
        void normalize(const Error<dim> &err);
        //friend class Display<dim>;
    };
*/

}

