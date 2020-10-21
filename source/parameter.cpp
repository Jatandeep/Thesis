#include"../include/parameter.h"
#include <deal.II/grid/grid_in.h>


using namespace dealii;
using namespace thesis;
using namespace parameters;

void GeometryModel::declare_param(ParameterHandler &prm){
    prm.enter_subsection("Geometry Model");
    {
        prm.declare_entry("Mesh File", "grid.inp", Patterns::Anything());
        prm.declare_entry("Global refinement", "0", Patterns::Integer());
        prm.declare_entry("Local refinement", "0", Patterns::Integer());
        prm.declare_entry("Grid scale", "1", Patterns::Double(0));
        prm.declare_entry("Plate dim", "1", Patterns::Double(0));
        prm.declare_entry("Crack length", "1", Patterns::Double(0));
        prm.declare_entry("Ref region height perc", "1", Patterns::Double(0));
        prm.declare_entry("Crack tip back ref perc", "1", Patterns::Double(0));
    }
    prm.leave_subsection();
}

void GeometryModel::parse_param(ParameterHandler &prm){
    prm.enter_subsection("Geometry Model");
    {
    meshfile = prm.get("Mesh File");
	gl_ref = prm.get_integer("Global refinement");
    lc_ref = prm.get_integer("Local refinement");
   	grid_scale = prm.get_double("Grid scale");
    b = prm.get_double("Plate dim");
    a = prm.get_double("Crack length");
    h = prm.get_double("Ref region height perc");
    x = prm.get_double("Crack tip back ref perc");
    }
    prm.leave_subsection();
}

void MaterialModel::declare_param(ParameterHandler &prm){
    prm.enter_subsection("Material Model");
    {
        prm.declare_entry("Lambda", "1", Patterns::Double(0));
        prm.declare_entry("Mu", "1", Patterns::Double(0));
    }
    prm.leave_subsection();
}

void MaterialModel::parse_param(ParameterHandler &prm){
    prm.enter_subsection("Material Model");
    {
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
        prm.declare_entry("Residual tolerance u", "1", Patterns::Double(0));
        prm.declare_entry("Residual tolerance d", "1", Patterns::Double(0));
        prm.declare_entry("Newton update tolerance u", "1", Patterns::Double(0));
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
        prm.declare_entry("Delta time initial", "0.1", Patterns::Double(0));
        prm.declare_entry("Delta time final", "0.1", Patterns::Double(0));
	    prm.declare_entry("Time change interval", "1", Patterns::Anything());
        prm.declare_entry("Time tolerance", "1", Patterns::Anything());
        prm.declare_entry("Output frequency", "5", Patterns::Integer());
        prm.declare_entry("Time adaptivity","true",Patterns::Selection("true|false"));
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
        delta_t = prm.get_double("Delta time initial");
   	    delta_t_f = prm.get_double("Delta time final");
	    time_change_interval = prm.get_double("Time change interval");
        time_tol = prm.get_double("Time tolerance");
        op_freq = prm.get_double("Output frequency");
        time_adap = prm.get("Time adaptivity");
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
        prm.declare_entry("Viscosity", "1", Patterns::Double(0));

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
        viscosity = prm.get_double("Viscosity");

    }
    prm.leave_subsection();
}

void TestCase::declare_param(ParameterHandler &prm){
    prm.enter_subsection("TestCase");
    {
        prm.declare_entry("Test case","tension",Patterns::Selection("tension|shear"));
    }
    prm.leave_subsection();
}
void TestCase::parse_param(ParameterHandler &prm){
    prm.enter_subsection("TestCase");
    {
        test = prm.get("Test case");
    }
    prm.leave_subsection();
}

void BoundaryConditions::declare_param(ParameterHandler &prm){
    prm.enter_subsection("BoundaryConditions");
    {
        prm.declare_entry("Tension x axis bottom","fixed",Patterns::Selection("fixed|free"));
        prm.declare_entry("Tension x axis top","free",Patterns::Selection("fixed|free"));
    }
    prm.leave_subsection();
}
void BoundaryConditions::parse_param(ParameterHandler &prm){
    prm.enter_subsection("BoundaryConditions");
    {
        uxb = prm.get("Tension x axis bottom");
        uxt = prm.get("Tension x axis top");
    }
    prm.leave_subsection();
}

void ModelingStrategy::declare_param(ParameterHandler &prm){
    prm.enter_subsection("ModelingStrategy");
    {
        prm.declare_entry("Initial crack strategy","M_I",Patterns::Selection("M_I|P_I|M_Id"));
        prm.declare_entry("Computing strategy","StandardNum",Patterns::Selection("StandardNum|lefm"));
        prm.declare_entry("Target factor fracture toughness", "2", Patterns::Double(0));
        prm.declare_entry("Target steps fracture toughness", "1000", Patterns::Double(0));
    }
    prm.leave_subsection();
}
void ModelingStrategy::parse_param(ParameterHandler &prm){
    prm.enter_subsection("ModelingStrategy");
    {
        strategy = prm.get("Initial crack strategy");
        comp_strategy = prm.get("Computing strategy");
        fac_ft = prm.get_double("Target factor fracture toughness");
        steps_ft = prm.get_double("Target steps fracture toughness");

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
    TestCase::declare_param(prm);
    BoundaryConditions::declare_param(prm);
    ModelingStrategy::declare_param(prm);
}

void AllParameters::parse_param(ParameterHandler &prm){
    geometrymodel.parse_param(prm);
    materialmodel.parse_param(prm);
    fesys.parse_param(prm);
    newtonraphson.parse_param(prm);
    linearsolver.parse_param(prm);
    time.parse_param(prm);
    pf.parse_param(prm);
    test_case.parse_param(prm);
    bc.parse_param(prm);
    mod_strategy.parse_param(prm);
}

