#include"../include/parameter.h"
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>
#include <deal.II/base/parameter_handler.h>


using namespace dealii;
using namespace thesis;
using namespace parameters;

void Geometrymodel::declare_param(ParameterHandler &prm){
    prm.enter_subsection("Geometry Model");
    {
        prm.declare_entry("Name", "gmsh", Patterns::Anything());
        prm.declare_entry("Mesh File", "grid.msh", Patterns::Anything());
    }
    prm.leave_subsection();
}

void Geometrymodel::parse_param(ParameterHandler &prm){
    prm.enter_subsection("Geometry Model");
    {
        name = prm.get("Name");
        meshfile = prm.get("Mesh File");
    }
    prm.leave_subsection();
}

void Materialmodel::declare_param(ParameterHandler &prm){
    prm.enter_subsection("Material Model");
    {
        prm.declare_entry("Elastic modulus", "1", Patterns::Anything());
        prm.declare_entry("Poisson ratio", "0.3", Patterns::Anything());
    }
    prm.leave_subsection();
}

void Materialmodel::parse_param(ParameterHandler &prm){
    prm.enter_subsection("Material Model");
    {
        elastic_mod = prm.get_double("Elastic modulus");
        poisson_ratio = prm.get_double("Poisson ratio");
    }
    prm.leave_subsection();
}


void FESys::declare_param(ParameterHandler &prm){
    prm.enter_subsection("Linear Elastic");
    {
        prm.declare_entry("Finite element degree", "1", Patterns::Integer());
        prm.declare_entry("Quad order", "1", Patterns::Integer(0));
        prm.declare_entry("Global refinement", "0", Patterns::Integer());
        prm.declare_entry("Max no of steps for solver", "1", Patterns::Integer(0));
        prm.declare_entry("CG Solver tolerance", "1", Patterns::Anything());
        prm.declare_entry("Relaxation parameter", "1", Patterns::Double(0));
        prm.declare_entry("Lambda", "1", Patterns::Double(0));
        prm.declare_entry("Mu", "1", Patterns::Double(0));
        prm.declare_entry("Act_ref", "0.3", Patterns::Double(0));
        prm.declare_entry("Act_cors", "0.03", Patterns::Double(0));
        prm.declare_entry("Starting time", "0", Patterns::Integer());
        prm.declare_entry("End time", "1", Patterns::Integer(0));
        prm.declare_entry("Delta time", "0.1", Patterns::Double(0));
        prm.declare_entry("Time tolerance", "1", Patterns::Anything());
        prm.declare_entry("Max newton iterations", "10", Patterns::Anything());
        prm.declare_entry("Residual tolerance", "1", Patterns::Double(0));



    }
    prm.leave_subsection();
}

void FESys::parse_param(ParameterHandler &prm){
    prm.enter_subsection("Linear Elastic");
    {
        fe_degree = prm.get_integer("Finite element degree");
        quad_order = prm.get_double("Quad order");
        gl_ref = prm.get_integer("Global refinement");
        cg_steps = prm.get_double("Max no of steps for solver");
        cg_tol = prm.get_double("CG Solver tolerance");
        relax_prm = prm.get_double("Relaxation parameter");
        lambda = prm.get_double("Lambda");
        mu = prm.get_double("Mu");
        act_ref = prm.get_double("Act_ref");
        act_cors = prm.get_double("Act_cors");
        start_time = prm.get_integer("Starting time");
        end_time = prm.get_integer("End time");
        delta_t = prm.get_double("Delta time");
        time_tol = prm.get_double("Time tolerance");
        max_new_ite = prm.get_integer("Max newton iterations");
        res_tol = prm.get_double("Residual tolerance");


    }
    prm.leave_subsection();
}

AllParameters::AllParameters(const std::string &filename){
    ParameterHandler prm;
    declare_param(prm);
    prm.parse_input(filename);
    parse_param(prm);
}

void AllParameters::declare_param(ParameterHandler &prm){
    Geometrymodel::declare_param(prm);
    Materialmodel::declare_param(prm);
    FESys::declare_param(prm);
}

void AllParameters::parse_param(ParameterHandler &prm){
    geometrymodel.parse_param(prm);
    materialmodel.parse_param(prm);
    fesys.parse_param(prm);
}

