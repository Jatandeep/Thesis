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
        prm.declare_entry("Local refinement", "0", Patterns::Integer());
        prm.declare_entry("Act_ref", "0.3", Patterns::Double(0));
        prm.declare_entry("Act_cors", "0.03", Patterns::Double(0));
        prm.declare_entry("Grid scale", "1", Patterns::Double(0));
    }
    prm.leave_subsection();
}

void GeometryModel::parse_param(ParameterHandler &prm){
    prm.enter_subsection("Geometry Model");
    {
        meshfile = prm.get("Mesh File");
	gl_ref = prm.get_integer("Global refinement");
       	lc_ref = prm.get_integer("Local refinement");
	act_ref = prm.get_double("Act_ref");
        act_cors = prm.get_double("Act_cors");
   	grid_scale = prm.get_double("Grid scale");
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
        prm.declare_entry("Viscosity", "1", Patterns::Double(0));

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
        viscosity = prm.get_double("Viscosity");
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
        prm.declare_entry("Residual tolerance u", "1", Patterns::Double(0));
        prm.declare_entry("Residual tolerance d", "1", Patterns::Double(0));
        prm.declare_entry("Newton update tolerance u", "1", Patterns::Double(0));
        prm.declare_entry("Newton update tolerance d", "1", Patterns::Double(0));
    }
    prm.leave_subsection();
}

void NewtonRaphson::parse_param(ParameterHandler &prm){
    prm.enter_subsection("Newton Raphson");
    {
        max_new_ite = prm.get_integer("Max newton iterations");
        res_tol_u = prm.get_double("Residual tolerance u");
        res_tol_d = prm.get_double("Residual tolerance d");
        nu_tol_u = prm.get_double("Newton update tolerance u");
        nu_tol_d = prm.get_double("Newton update tolerance d");
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
        prm.declare_entry("Starting time", "0", Patterns::Double(0));
        prm.declare_entry("End time", "1", Patterns::Double(0));
        prm.declare_entry("Delta time", "0.1", Patterns::Double(0));
        prm.declare_entry("Delta time final", "0.1", Patterns::Double(0));
	prm.declare_entry("Time change interval", "1", Patterns::Anything());
        prm.declare_entry("Time tolerance", "1", Patterns::Anything());
        prm.declare_entry("Output frequency", "5", Patterns::Integer());
        prm.declare_entry("Alpha", "1", Patterns::Double(0));
        prm.declare_entry("Beta", "1", Patterns::Double(0));
    }
    prm.leave_subsection();
}

void Time::parse_param(ParameterHandler &prm){
    prm.enter_subsection("Time");
    {
        start_time = prm.get_double("Starting time");
        end_time = prm.get_double("End time");
        delta_t = prm.get_double("Delta time");
   	delta_t_f = prm.get_double("Delta time final");
	time_change_interval = prm.get_double("Time change interval");
        time_tol = prm.get_double("Time tolerance");
        op_freq = prm.get_double("Output frequency");
        alpha = prm.get_double("Alpha");
        beta = prm.get_double("Beta");
    }
    prm.leave_subsection();
}


void PhaseFieldMethod::declare_param(ParameterHandler &prm){
    prm.enter_subsection("PhaseField");
    {
        prm.declare_entry("Critical energy release rate", "1", Patterns::Double(0));
        prm.declare_entry("Length scale parameter", "1", Patterns::Double(0));
	prm.declare_entry("Small positive parameter", "1", Patterns::Double(0));
	prm.declare_entry("Total displacement", "1", Patterns::Double(0));
    }
    prm.leave_subsection();
}

void PhaseFieldMethod::parse_param(ParameterHandler &prm){
    prm.enter_subsection("PhaseField");
    {
        g_c = prm.get_double("Critical energy release rate");
        l = prm.get_double("Length scale parameter");
        k = prm.get_double("Small positive parameter");
	u_total = prm.get_double("Total displacement");
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
    PhaseFieldMethod::declare_param(prm);
}

void AllParameters::parse_param(ParameterHandler &prm){
    geometrymodel.parse_param(prm);
    materialmodel.parse_param(prm);
    fesys.parse_param(prm);
    newtonraphson.parse_param(prm);
    linearsolver.parse_param(prm);
    time.parse_param(prm);
    pf.parse_param(prm);
}

