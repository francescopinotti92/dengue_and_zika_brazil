//
//  main.cpp
//  CrossFlavivirus
//
//  Created by MacBook Pro on 05/01/2021.
//

#include <iostream>
#include <string>
#include <nlohmann/json.hpp>
#include "event_vec.hpp"
#include "individual.hpp"
#include "utils.hpp"
#include <chrono>


int main(int argc, const char * argv[]) {

    // arguments from bash
    
    const char* pathParams = argv[1];
    int startingSim = atoi( argv[2] );
    int Nsim = atoi( argv[3] );
    int seed = atoi( argv[4] );
    
    // read model parameters
    
    //string path_to_params = "/Users/user/Documents/project_simulations/CrossFlavivirus/test1/parameters/settings_0.json";
    //std::ifstream readFile(path_to_params);
    
    std::ifstream readFile( pathParams );
    nlohmann::json params = nlohmann::json::parse(readFile);
    readFile.close();
    
    int scenario_id = params["scenario_id"];                        // ID scenario
    int N = params["N"];                                            // Pop size
    int Tintro = params["Tintro"];                                  // Timing Zika introduction
    int Tsample = params["Tsample"];                                // Timing begin sampling
    int Tmax = params["Tmax"];
    int nstrains = params["nstrains"];
    vector<int> whichDs = params["whichDs"].get<vector<int>>();
    double intro_rate = params["intro_rate"];
    double mean_degree = params["mean_degree"];
    double seasonal_force_D = params["seasonal_force_D"];
    double seasonal_force_Z = params["seasonal_force_Z"];
    double seasonal_max_t_D = params["seasonal_max_t_D"];
    double seasonal_max_t_Z = params["seasonal_max_t_Z"];
    double age_shape = params["age_shape"];
    double age_max_prob = params["age_max_prob"];
    vector<double> trs = params["transmissibility"];
    vector<double> sigs = params["sigmas"];
    vector<double> ltc_probs = params["epsilons"];

    
    //double rho11 = params["rho11"];
    //double rho12 = params["rho12"];
    //double rho21 = params["rho21"];
    
    double lD = static_cast<double>(params["lD"]) * 365.; // assume this is given in years => convert to days
    double lZ = static_cast<double>(params["lZ"]) * 365.; // assume this is given in years => convert to days
    double rhoDZ = params["rhoDZ"];
    double rhoZD = params["rhoZD"];
    double rhoZD0 = params["rhoZD0"];
    double csiDZ = params["csiDZ"];
    double csiZD = params["csiZD"];
    double gamma1D = params["gamma1D"];
    double gamma1Z = params["gamma1Z"];
    double gamma2D = params["gamma2D"];
    double gamma2Z = params["gamma2Z"];
    
    vector<int> endemic_strains = {};
    vector<int> all_strains = {};
    
    // create folders if they do not exist
    
    string save_path = params["save_path"];
    if (exists_path(save_path)==0) {
        mkdir(save_path.c_str(), ACCESSPERMS);
    }
    
    ofstream write_prev;
    ofstream write_detailed_events;
    ofstream write_AR_age_endemic;
    ofstream write_AR_age_invasive;
    //ofstream write_age;
    ofstream write_incidence;
    ofstream write_seroprevalenceD;
    
    // output flags
    
    //bool writeAge = params["outputFlags"]["writeAge"].get<bool>();
    bool writePrev = params["outputFlags"]["writePrevalence"].get<bool>();
    bool writeInc = params["outputFlags"]["writeIncidence"].get<bool>();
    bool writeAR_D = params["outputFlags"]["writeAREndemic"].get<bool>();
    bool writeAR_Z = params["outputFlags"]["writeARInvasive"].get<bool>();
    bool writeSeroPrevD = params["outputFlags"]["writeSeroPrevalenceEndemic"].get<bool>();
    bool writeInfEvs = params["outputFlags"]["writeInfectionEvents"].get<bool>();
    
    
    mt19937_64 mymt = mt19937_64(1000);
    mymt.seed( seed );
    assign_mt(&mymt);
    
    CompartmentalModel compM = CompartmentalModel( nstrains, trs, sigs, ltc_probs );
    compM.set_immunity_params( lD, lZ, rhoDZ, rhoZD, rhoZD0, csiDZ, csiZD, gamma1D, gamma1Z, gamma2D, gamma2Z );
    
    for ( int i = 0; i < nstrains; i++ ) {
        if ( compM.is_endemic( i ) )
            endemic_strains.push_back( i );
        all_strains.push_back( i );
    }
    
    std::chrono::steady_clock::time_point begint = std::chrono::steady_clock::now();

    for ( int sim = startingSim; sim < startingSim + Nsim; ++sim ) {
        compM.clear_incidence();
        Community comm = Community(N, intro_rate, mean_degree, seasonal_force_D, seasonal_force_Z, seasonal_max_t_D, seasonal_max_t_Z, age_shape, age_max_prob);
        
        if ( writePrev )
            write_prev.open(save_path + "/prevalence_" + to_string(scenario_id) + "_" + to_string(sim) + ".txt");
        
        if ( writeAR_D )
            write_AR_age_endemic.open(save_path + "/AR_age_endemic_" + to_string(scenario_id) + "_" + to_string(sim) + ".txt");
        
        if ( writeAR_Z )
            write_AR_age_invasive.open(save_path + "/AR_age_invasive_" + to_string(scenario_id) + "_" + to_string(sim) + ".txt");
        
        if ( writeInc )
            write_incidence.open(save_path + "/incidence_" + to_string(scenario_id) + "_" + to_string(sim) + ".txt");

        if ( writeSeroPrevD )
            write_seroprevalenceD.open(save_path + "/seroprevalenceD_" + to_string(scenario_id) + "_" + to_string(sim) + ".txt");

        //if ( writeAge )
        //    string write_age_path = save_path + "/age_" + to_string(scenario_id) + "_" + to_string(sim) + ".txt";
        
        string write_tlastD_path = save_path + "/tLastDengueInfection_" + to_string(scenario_id) + "_" + to_string(sim) + ".txt";

        if ( writeInfEvs ) {
            write_detailed_events.open(save_path + "/all_evs_" + to_string(scenario_id) + "_" + to_string(sim) + ".txt");
            write_detailed_events << "t,strain,age,nD,nZ,tD,t1D,tZ" << std::endl;
        }
        
        int t = 0;
        
        // initialize high levels of immunity to endemic strains in population
        for ( int strain: whichDs )
            comm.initialize_immunity_strain(strain, compM);

        for (; t < Tmax; ++t) {
                        
            if (t == Tintro) {
                //comm.print_TlastDengueInfection(t, compM, write_tlastD_path);
                comm.set_random_infection( compM.get_Z_id(), 10, t, compM);            // introduce Zika
                //comm.print_age_distribution(t, write_age_path); // write age distribution
                if ( Tintro == Tsample )
                    compM.clear_incidence(); // clear incidence now (will save space on hard disk slightly)
            }
            
            if ( ( t % 20 == 0 ) and ( writeInc ) ) {
                compM.print_cumulative_incidence(t, write_incidence);
            }
            
            if ( ( t % 20 == 0 ) and ( writePrev ) ) {
                vector<int> hist = comm.get_prevalence_strains(compM);
                write_prev << t << "," << stringify_vector(hist) << std::endl;

            }
            
            if ( t < Tsample )
                comm.do_dynamics_step( t, compM );                        // epi dynamics
            else {
                comm.do_dynamics_step_output( t, compM );                 // epi dynamics with recording
                
                if ( ( t > Tsample ) and ( t % 365 == 0 ) ) {
                    if ( writeAR_D )
                        comm.print_AR_age_endemic_expZ(t, compM, write_AR_age_endemic);
                    if ( writeAR_Z )
                        comm.print_AR_age_invasive_exp(t, compM, write_AR_age_invasive);
                    if ( writeInfEvs )
                        comm.print_infection_events_flush( t, compM, write_detailed_events );
                }
                
                /*
                if ( ( t == Tmax - 1 ) & ( writeInfEvs ) ) {
                    comm.print_infection_events( t, compM, write_all_inf_evs_path );
                    comm.flush_events();
                }*/

            
                if ( ( t >= Tsample ) and ( t % (10 * 365) == 0 ) and ( writeSeroPrevD ) ) {
                    comm.print_seroprevalenceD(t, compM, write_seroprevalenceD);
                }
                
            }
            
            if ( t < Tintro )
                comm.do_demography_step(t, compM, false, false);
            else {
                if ( t < Tsample )
                    comm.do_demography_step(t, compM, false, false);
                else
                    comm.do_demography_step(t, compM, false, true);
            }

            comm.do_intro_step(t, whichDs, compM);          // intro step
        }
        
        if ( writePrev )
            write_prev.close();
        
        if ( writeAR_D )
            write_AR_age_endemic.close();
        
        if ( writeAR_Z )
            write_AR_age_invasive.close();

        //write_detailed_events.close();
        
        if ( writeInc )
            write_incidence.close();
        
        if ( writeSeroPrevD )
            write_seroprevalenceD.close();
        
        if ( writeInfEvs )
            write_detailed_events.close();
        
        comm.flush_events();
        compM.clear_incidence();
                
    }

    std::chrono::steady_clock::time_point endt = std::chrono::steady_clock::now();
    std::cout << "Time elapsed = " << std::chrono::duration_cast<std::chrono::microseconds>(endt - begint).count()/1000000 << "[s]" << std::endl;

}
