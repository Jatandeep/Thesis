#ifndef HYPERCUBEWITHREFINEDHOLE_H
#define HYPERCUBEWITHREFINEDHOLE_H

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>
#include <iostream>
#include <fstream>
#include <cmath>


using namespace dealii;

namespace HyperCubeWithRefinedHole
{

//This function refines the mesh around the inner region a number of times 
//according to the given number in the input
template<int dim>
void set_and_execute_refinements(unsigned int number_refinements, Triangulation<dim> &triangulation, types::boundary_id boundary_id_inner_hole)
{
	for(unsigned int i=0; i<number_refinements; i++)
	{
		typename Triangulation<dim>::active_cell_iterator cell= triangulation.begin_active(),
												endc = triangulation.end();
		for(; cell!=endc; ++cell)
		{
			for(unsigned int j=0; j<GeometryInfo<dim>::faces_per_cell; j++)
			{
				if (cell->face(j)->at_boundary() && cell->face(j)->boundary_id() == boundary_id_inner_hole)
				{
					cell->set_refine_flag();
					break;
				}
			}
		}
		triangulation.execute_coarsening_and_refinement();		
	}
}



template<int dim>
void set_boundary_ids(Triangulation<dim> &triangulation, unsigned int id_Dirichlet, unsigned int id_Neumann, double outer_radius)
{
	typename Triangulation<dim>::active_cell_iterator cell= triangulation.begin_active(),
											endc = triangulation.end();

	for(; cell!=endc; ++cell)
	{
		for(unsigned int j=0; j<GeometryInfo<dim>::faces_per_cell; j++)
		{
			if (cell->face(j)->at_boundary() && (std::abs((cell->face(j)->center()[0] - outer_radius)) < 1e-10) )
			{
				cell->face(j)->set_boundary_id(id_Neumann);
			}
			else if (cell->face(j)->at_boundary() && std::abs((cell->face(j)->center()[0] + outer_radius)) < 1e-10)
			{
				cell->face(j)->set_boundary_id(id_Dirichlet);
			}
		}
		
	}	
}
//This function uses the GridGenerator to generate a
//hypercupe with zylindrical hole
template<int dim>
void generate_grid(Triangulation<dim> &triangulation,
							unsigned int int_nbr_refinements,
							unsigned int int_nbr_global_refinements,
							unsigned int id_dirichelt_boundary,
						  unsigned int id_neumann_boundary )
{
	const double outer_radius = 1.0;
	const double inner_radius = 0.5;
	const Point<dim> center;

	GridGenerator::hyper_cube_with_cylindrical_hole(triangulation,
													inner_radius,
												 outer_radius,
												 0.5,
												 1,
												 false /*boundary_id_inner_hole is set to 1*/);
	
	//create both manifold descriptions for the two and three dimensional case
	const SphericalManifold<dim> manifold_description_2d(center);
	//const CylindricalManifold<dim> manifold_description_3d(dim-1);
	types::boundary_id boundary_id_inner_hole=1;
	types::manifold_id manifold_id_inner_hole=1;
	//Set the manifold id of all boundary faces and edges with given boundary id 
	triangulation.set_all_manifold_ids_on_boundary(boundary_id_inner_hole,manifold_id_inner_hole);
	/*If this is not done and the manifold_id equals number::invalid_manifold_id (which is default)
	*the triangulation object queries the boundary_id if the face is at the boundary or the material_id
	*/
	//Checker wether or not the mesh is for dim==2 or dim==3
	if(dim==2)
	{
		triangulation.set_manifold (manifold_id_inner_hole, manifold_description_2d);
	}/*
	else if(dim==3)
	{
		triangulation.set_manifold (manifold_id_inner_hole, manifold_description_3d);
	}
	*/
	triangulation.refine_global(int_nbr_global_refinements);
	set_and_execute_refinements(int_nbr_refinements, triangulation, boundary_id_inner_hole);
	triangulation.set_manifold(manifold_id_inner_hole/*reset to default FlatManifold*/);
		
	set_boundary_ids(triangulation,id_dirichelt_boundary, id_neumann_boundary, outer_radius);


}




}//end of namespace
#endif
