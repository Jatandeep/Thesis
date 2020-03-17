#include"../include/parameter.h"
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>
#include <deal.II/base/parameter_handler.h>


using namespace dealii;
using namespace thesis;
using namespace parameters;

void GeometryModel::declare_param(ParameterHandler &prm){
    prm.enter_subsection("Geometry Model");
    {
        prm.declare_entry("Mesh File", "grid.inp", Patterns::Anything());
        prm.declare_entry("Global refinement", "0", Patterns::Integer());
        prm.declare_entry("Act_ref", "0.3", Patterns::Double(0));
        prm.declare_entry("Act_cors", "0.03", Patterns::Double(0));
    }
    prm.leave_subsection();
}

void GeometryModel::parse_param(ParameterHandler &prm){
    prm.enter_subsection("Geometry Model");
    {
        meshfile = prm.get("Mesh File");
	gl_ref = prm.get_integer("Global refinement");
        act_ref = prm.get_double("Act_ref");
        act_cors = prm.get_double("Act_cors");
    }
    prm.leave_subsection();
}

void MaterialModel::declare_param(ParameterHandler &prm){
    prm.enter_subsection("Material Model");
    {
        prm.declare_entry("Elastic modulus", "1", Patterns::Anything());
        prm.declare_entry("Poisson ratio", "0.3", Patterns::Anything());
        prm.declare_entry("Lambda", "1", Patterns::Double(0));
        prm.declare_entry("Mu", "1", Patterns::Double(0));
    }
    prm.leave_subsection();
}

void MaterialModel::parse_param(ParameterHandler &prm){
    prm.enter_subsection("Material Model");
    {
        elastic_mod = prm.get_double("Elastic modulus");
        poisson_ratio = prm.get_double("Poisson ratio");
        lambda = prm.get_double("Lambda");
        mu = prm.get_double("Mu");
    }
    prm.leave_subsection();
}


void FESys::declare_param(ParameterHandler &prm){
    prm.enter_subsection("FE System");
    {
        prm.declare_entry("Finite element degree", "1", Patterns::Integer());
        prm.declare_entry("Quad order", "1", Patterns::Integer(0));
    }
    prm.leave_subsection();
}

void FESys::parse_param(ParameterHandler &prm){
    prm.enter_subsection("FE System");
    {
        fe_degree = prm.get_integer("Finite element degree");
        quad_order = prm.get_double("Quad order");
    }
    prm.leave_subsection();
}


void NewtonRaphson::declare_param(ParameterHandler &prm){
    prm.enter_subsection("Newton Raphson");
    {
        prm.declare_entry("Max newton iterations", "10", Patterns::Anything());
        prm.declare_entry("Residual tolerance", "1", Patterns::Double(0));
    }
    prm.leave_subsection();
}

void NewtonRaphson::parse_param(ParameterHandler &prm){
    prm.enter_subsection("Newton Raphson");
    {
        max_new_ite = prm.get_integer("Max newton iterations");
        res_tol = prm.get_double("Residual tolerance");
    }
    prm.leave_subsection();
}

void LinearSolver::declare_param(ParameterHandler &prm){
    prm.enter_subsection("Linear Solver");
    {
        prm.declare_entry("CG Solver tolerance", "1", Patterns::Anything());
        prm.declare_entry("Relaxation parameter", "1", Patterns::Double(0));
    }
    prm.leave_subsection();
}

void LinearSolver::parse_param(ParameterHandler &prm){
    prm.enter_subsection("Linear Solver");
    {
        cg_tol = prm.get_double("CG Solver tolerance");
        relax_prm = prm.get_double("Relaxation parameter");
    }
    prm.leave_subsection();
}


void Time::declare_param(ParameterHandler &prm){
    prm.enter_subsection("Time");
    {
        prm.declare_entry("Starting time", "0", Patterns::Integer());
        prm.declare_entry("End time", "1", Patterns::Integer(0));
        prm.declare_entry("Delta time", "0.1", Patterns::Double(0));
        prm.declare_entry("Time tolerance", "1", Patterns::Anything());
    }
    prm.leave_subsection();
}

void Time::parse_param(ParameterHandler &prm){
    prm.enter_subsection("Time");
    {
        start_time = prm.get_integer("Starting time");
        end_time = prm.get_integer("End time");
        delta_t = prm.get_double("Delta time");
        time_tol = prm.get_double("Time tolerance");
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
    GeometryModel::declare_param(prm);
    MaterialModel::declare_param(prm);
    FESys::declare_param(prm);
    NewtonRaphson::declare_param(prm);
    LinearSolver::declare_param(prm);
    Time::declare_param(prm);
}

void AllParameters::parse_param(ParameterHandler &prm){
    geometrymodel.parse_param(prm);
    materialmodel.parse_param(prm);
    fesys.parse_param(prm);
	newtonraphson.parse_param(prm);
	linearsolver.parse_param(prm);
	time.parse_param(prm);
}

