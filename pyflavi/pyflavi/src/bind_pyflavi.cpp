//
//  bind_pyFlavi.cpp
//  CrossFlavivirus
//
//  Created by user on 24/08/2021.
//

#include <stdio.h>

#include "pyflavi.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

// module: is the name of the python module we are creating, should coincide with cpp file name?
// python complained when I used 'mymodule' instead of 'module'
// m: is an instance of py::module_, necessary to create bindings

namespace py = pybind11;

PYBIND11_MODULE(pyflavi, m) {
    m.doc() = "python binding for c++ code simulating flavivirus dynamics"; // optional module docstring

    //py::bind_vector<std::vector<int>>(m, "IntVector3D");
    //py::bind_vector<std::vector<std::vector<std::vector<std::vector<int>>>>>(m, "IntVector4D");
    py::bind_vector<std::vector<std::vector<std::vector<std::vector<int>>>>>(m, "IntVector4D");
    
    py::class_<IndividualInfectionEvent>( m, "IndividualInfectionEvent" ).def_readwrite( "strain", &IndividualInfectionEvent::strain).def_readwrite( "t", &IndividualInfectionEvent::t );
    
    py::bind_vector<std::vector<IndividualInfectionEvent>>(m, "IndividualInfectionEventVec");

    
    py::class_<EventInfectionO>( m, "EventInfectionO").def_readwrite("id", &EventInfectionO::id).def_readwrite("t", &EventInfectionO::t).def_readwrite("age", &EventInfectionO::age).def_readwrite("strain", &EventInfectionO::strain).def_readwrite("nD", &EventInfectionO::previousExpD ).def_readwrite("nZ", &EventInfectionO::previousExpZ ).def_readwrite("past_events", &EventInfectionO::past_evs );
    
    py::bind_vector<std::vector<std::vector<EventInfectionO>>>(m, "EventInfectionOVec");

    m.def("simulate_DENV", &simulate_DENV, "Simulates spread",
          py::arg("Nsim"),
          py::arg("seed"),
          py::arg("Tmax"),
          py::arg("Tsample"),
          py::arg("N"),
          py::arg("nstrains"),
          py::arg("beta"),
          py::arg("epsi"),
          py::arg("sigma"),
          py::arg("mean_degree"),
          py::arg("lD"),
          py::arg("gammaDD1"),
          py::arg("gammaDD2"),
          py::arg("intro_rate"),
          py::arg("seasonal_force"),
          py::arg("seasonal_peak"),
          py::arg("age_shape"),
          py::arg("age_max_prob"),
          py::arg("whichDs"),
          py::arg("doImmunityIntro"),
          py::arg("doDemography"),
          py::arg("checkpts_incidence"),
          py::arg("checkpts_seroprevalence") );
    
    m.def("simulate_DENV_detailed_incidence", &simulate_DENV_detailed_incidence, "Simulates spread",
          py::arg("Nsim"),
          py::arg("seed"),
          py::arg("Tmax"),
          py::arg("Tsample"),
          py::arg("N"),
          py::arg("nstrains"),
          py::arg("beta"),
          py::arg("epsi"),
          py::arg("sigma"),
          py::arg("mean_degree"),
          py::arg("lD"),
          py::arg("gammaDD1"),
          py::arg("gammaDD2"),
          py::arg("intro_rate"),
          py::arg("seasonal_force"),
          py::arg("seasonal_peak"),
          py::arg("age_shape"),
          py::arg("age_scale"),
          py::arg("whichDs"),
          py::arg("checkpts") );
    
    m.def("simulate_DENV_sero", &simulate_DENV_sero, "Simulates spread",
          py::arg("Nsim"),
          py::arg("seed"),
          py::arg("Tmax"),
          py::arg("Tsample"),
          py::arg("N"),
          py::arg("nstrains"),
          py::arg("beta"),
          py::arg("epsi"),
          py::arg("sigma"),
          py::arg("mean_degree"),
          py::arg("lD"),
          py::arg("gammaDD1"),
          py::arg("gammaDD2"),
          py::arg("intro_rate"),
          py::arg("seasonal_force"),
          py::arg("seasonal_peak"),
          py::arg("age_shape"),
          py::arg("age_scale"),
          py::arg("sero_lam1"),
          py::arg("sero_lam2"),
          py::arg("sero_a0"),
          py::arg("checkpts_incidence"),
          py::arg("checkpts_seroprevalence"));
    
    m.def("simulate_ZIKV_sero_detailed", &simulate_ZIKV_sero_detailed, "Simulates spread",
          py::arg("Nsim"),
          py::arg("seed"),
          py::arg("Tmax"),
          py::arg("TintroZ"),
          py::arg("N"),
          py::arg("nstrains"),
          py::arg("betaD"),
          py::arg("betaZ"),
          py::arg("epsiD"),
          py::arg("epsiZ"),
          py::arg("sigmaD"),
          py::arg("sigmaZ"),
          py::arg("mean_degree"),
          py::arg("lD"),
          py::arg("lZ"),
          py::arg("intro_rate"),
          py::arg("seasonal_forceD"),
          py::arg("seasonal_forceZ"),
          py::arg("seasonal_peakD"),
          py::arg("seasonal_peakZ"),
          py::arg("rhoDZ"),
          py::arg("rhoZD"),
          py::arg("rhoZD0"),
          py::arg("csiDZ"),
          py::arg("csiZD"),
          py::arg("gamma1D"),
          py::arg("gamma1Z"),
          py::arg("gamma2D"),
          py::arg("gamma2Z"),
          py::arg("age_shape"),
          py::arg("age_scale"),
          py::arg("sero_lam1"),
          py::arg("sero_lam2"),
          py::arg("sero_a0") );
    
    
    
    m.def("simulate_DENV_removal", &simulate_DENV_removal, "Simulates spread",
          py::arg("Nsim"),
          py::arg("seed"),
          py::arg("Tmax"),
          py::arg("Tremoval"),
          py::arg("dT"),
          py::arg("N"),
          py::arg("nstrains"),
          py::arg("beta"),
          py::arg("epsi"),
          py::arg("sigma"),
          py::arg("mean_degree"),
          py::arg("lD"),
          py::arg("intro_rate"),
          py::arg("seasonal_force"),
          py::arg("seasonal_peak"),
          py::arg("age_shape"),
          py::arg("age_max_prob"),
          py::arg("whichDs"),
          py::arg("checkpts"));

    m.def("simulate_ZIKV", &simulate_ZIKV, "Simulates both ZIKV and DENV",
          py::arg("Nsim"),
          py::arg("seed"),
          py::arg("Tmax"),
          py::arg("Tsample"),
          py::arg("TintroZ"),
          py::arg("N"),
          py::arg("nstrains"),
          py::arg("whichDs"),
          py::arg("checkpts"),
          py::arg("printZonly"),
          py::arg("betaD"),
          py::arg("betaZ"),
          py::arg("epsiD"),
          py::arg("epsiZ"),
          py::arg("sigmaD"),
          py::arg("sigmaZ"),
          py::arg("mean_degree"),
          py::arg("lD"),
          py::arg("lZ"),
          py::arg("intro_rate"),
          py::arg("seasonal_forceD"),
          py::arg("seasonal_forceZ"),
          py::arg("seasonal_peakD"),
          py::arg("seasonal_peakZ"),
          py::arg("age_shape"),
          py::arg("age_max_prob"),
          py::arg("rhoDZ"),
          py::arg("rhoZD"),
          py::arg("rhoZD0"),
          py::arg("csiDZ"),
          py::arg("csiZD"),
          py::arg("gamma1D"),
          py::arg("gamma1Z"),
          py::arg("gamma2D"),
          py::arg("gamma2Z"),
          py::arg("nZIKV0") = 0);
    
    m.def("simulate_ZIKV_sero", &simulate_ZIKV_sero, "Simulates both ZIKV and DENV",
          py::arg("Nsim"),
          py::arg("seed"),
          py::arg("Tmax"),
          py::arg("Tsample"),
          py::arg("TintroZ"),
          py::arg("N"),
          py::arg("nstrains"),
          py::arg("whichDs"),
          py::arg("checkpts"),
          py::arg("printZonly"),
          py::arg("betaD"),
          py::arg("betaZ"),
          py::arg("epsiD"),
          py::arg("epsiZ"),
          py::arg("sigmaD"),
          py::arg("sigmaZ"),
          py::arg("mean_degree"),
          py::arg("lD"),
          py::arg("lZ"),
          py::arg("intro_rate"),
          py::arg("seasonal_forceD"),
          py::arg("seasonal_forceZ"),
          py::arg("seasonal_peakD"),
          py::arg("seasonal_peakZ"),
          py::arg("age_shape"),
          py::arg("age_max_prob"),
          py::arg("rhoDZ"),
          py::arg("rhoZD"),
          py::arg("rhoZD0"),
          py::arg("csiDZ"),
          py::arg("csiZD"),
          py::arg("gamma1D"),
          py::arg("gamma1Z"),
          py::arg("gamma2D"),
          py::arg("gamma2Z"),
          py::arg("sero_lam1"),
          py::arg("sero_lam2"),
          py::arg("sero_a0"),
          py::arg("nZIKV0") = 0);
   
    m.def("simulate_ZIKV_parityDependentCrossReactivity", &simulate_ZIKV_parityDependentCrossReactivity, "blabla",
          py::arg("Nsim"),
          py::arg("seed"),
          py::arg("Tmax"),
          py::arg("Tsample"),
          py::arg("TintroZ"),
          py::arg("N"),
          py::arg("nstrains"),
          py::arg("whichDs"),
          py::arg("checkpts"),
          py::arg("printZonly"),
          py::arg("betaD"),
          py::arg("betaZ"),
          py::arg("epsiD"),
          py::arg("epsiZ"),
          py::arg("sigmaD"),
          py::arg("sigmaZ"),
          py::arg("mean_degree"),
          py::arg("lD"),
          py::arg("lZ"),
          py::arg("intro_rate"),
          py::arg("seasonal_forceD"),
          py::arg("seasonal_forceZ"),
          py::arg("seasonal_peakD"),
          py::arg("seasonal_peakZ"),
          py::arg("age_shape"),
          py::arg("age_max_prob"),
          py::arg("rhoZD0"),
          py::arg("rhoZD1"),
          py::arg("rhoZD2"),
          py::arg("chiZD0"),
          py::arg("chiZD1"),
          py::arg("chiZD2"));
    
    // use py::arg("name_arg") to enable named arguments in python
    // use py::arg("name_arg") = x to enable named arguments with default parameters
    /*
    m.def("add", &add, "A function which adds two numbers",
    py::arg("i") = 1, py::arg("j") = 2);
     */
}
