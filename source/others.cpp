#include "../include/ElasticProblem.h"
#include "../include/others.h"

using namespace dealii;
using namespace thesis;








template <int dim>
void Others<dim>::right_hand_side (const std::vector<Point<dim> > &points,
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

template class thesis::Others<2>;

//template class thesis::right_hand_side(const std::vector<dealii::Point<2> > &points,
//std::vector<dealii::Tensor<1, 2> >   &values);


//template class thesis::ElasticProblem<2>;
