#include "../include/ElasticProblem.h"
#include "../include/others.h"

using namespace dealii;
using namespace thesis;


template <int dim>
double Others<dim>::value(const Point<dim> &point, unsigned int c)const{}


template <int dim>
void Others<dim>::vector_value (const Point<dim> &points,
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

template class thesis::Others<2>;
template class thesis::Others<3>;
