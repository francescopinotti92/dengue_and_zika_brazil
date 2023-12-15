//
//  pyFlavi.cpp
//  pyFlavi
//
//  Created by user on 23/08/2021.
//

#include "pyflavi.hpp"


std::pair<std::vector<std::vector<std::vector<std::vector<int>>>>, std::vector<std::vector<std::vector<std::vector<int>>>>> simulate_DENV(const int& Nsim,
                   const int& seed,
                   const int& Tmax,
                   const int& Tsample,
                   const int& N,
                   const int& nstrains,
                   const double& beta,
                   const double& epsi,
                   const double& sigma,
                   const double& mean_degree,
                   const double& lD,
                   const double& gammaDD1,
                   const double& gammaDD2,
                   const double& intro_rate,
                   const double& seasonal_force,
                   const double& seasonal_peak,
                   const double& age_shape,
                   const double& age_max_prob,
                   const std::vector<int>& whichDs,
                   const bool& doImmunityIntro,
                   const bool& doDemography,
                   const std::vector<int>& checkpts_inc, const std::vector<int>& checkpts_sero) {
    
    std::vector<double> trs = std::vector<double>( nstrains, beta );
    std::vector<double> epsis = std::vector<double>( nstrains, epsi );
    std::vector<double> sigs = std::vector<double>( nstrains, sigma );
    
    double lZ = 365.;
    
    double rhoDZ = 1.;
    double rhoZD = 1.;
    double rhoZD0 = 1.;
    double csiDZ = 1.;
    double csiZD = 1.;
    double gamma1D = gammaDD1;
    double gamma1Z = 1.;
    double gamma2D = gammaDD2;
    double gamma2Z = 1.;
    
    std::vector<int> endemic_strains = {};
    std::vector<int> all_strains = {};
    std::vector<int> intro_strains = {};
        
    // output flags
    
    std::mt19937_64 mymt = std::mt19937_64(1000);
    mymt.seed( seed );
    assign_mt(&mymt);
    
    CompartmentalModel compM = CompartmentalModel( nstrains, trs, sigs, epsis );
    compM.set_immunity_params( lD * 365., lZ, rhoDZ, rhoZD, rhoZD0, csiDZ, csiZD, gamma1D, gamma1Z, gamma2D, gamma2Z );
    
    for ( int i = 0; i < nstrains; i++ ) {
        if ( compM.is_endemic( i ) )
            endemic_strains.push_back( i );
        all_strains.push_back( i );
    }
    
    for ( int strain: whichDs )
        intro_strains.push_back( strain );
   
    // prepare output
    int n_checkpts_inc = static_cast<int>( checkpts_inc.size() );
    assert( n_checkpts_inc > 0 );
    std::vector<std::vector<std::vector<std::vector<int>>>> pooled_histo;
    pooled_histo = std::vector<std::vector<std::vector<std::vector<int>>>>( Nsim, std::vector<std::vector<std::vector<int>>>( n_checkpts_inc, std::vector<std::vector<int>>( 4, std::vector<int>( 101, 0 ) ) ) );
    
    int n_checkpts_sero = static_cast<int>( checkpts_sero.size() );
    assert( n_checkpts_sero > 0 );
    std::vector<std::vector<std::vector<std::vector<int>>>> seroprevalence;
    seroprevalence = std::vector<std::vector<std::vector<std::vector<int>>>>( Nsim, std::vector<std::vector<std::vector<int>>>( n_checkpts_sero, std::vector<std::vector<int>>( 5, std::vector<int>( 101, 0 ) ) ) );
    
    
    for ( int sim = 0; sim < Nsim; ++sim ) {
        
       
        int id_next_checkpt_inc  = 0;
        int id_next_checkpt_sero = 0;
        
        int t_next_checkpt_inc   = checkpts_inc[ id_next_checkpt_inc   ];
        int t_next_checkpt_sero  = checkpts_sero[ id_next_checkpt_sero ];

            
        Community comm = Community(N,
                                   intro_rate,
                                   mean_degree,
                                   seasonal_force,
                                   0.,
                                   seasonal_peak,
                                   0 );
        
        comm.initialize_population_ageonly( age_max_prob, age_shape );

        
        // initialize high levels of immunity to endemic strains in population
        if ( doImmunityIntro ) {
            for ( int strain: intro_strains )
                comm.initialize_immunity_strain( strain, compM );
        }

        for ( int t = 0; t < Tmax; ++t ) {
                        
            if ( t == Tsample ) {
                compM.clear_incidence(); // clear incidence now (will save space on hard disk slightly)
            }
            
            if ( t < Tsample )
                comm.do_dynamics_step( t, compM );                        // epi dynamics
            else {
                comm.do_dynamics_step_output( t, compM );
                
                // output: seroprevalence
                if ( t == t_next_checkpt_sero ) {
                    
                    comm.fill_seroprevalenceD( t, seroprevalence[ sim ][ id_next_checkpt_sero ], compM );
                    
                    ++id_next_checkpt_sero;
    
                    if ( id_next_checkpt_sero < n_checkpts_sero )
                        t_next_checkpt_sero = checkpts_sero[id_next_checkpt_sero];
                    else
                        t_next_checkpt_sero = Tmax + 1;
                    
                }
                
                
                // output: incidence
                if ( t == t_next_checkpt_inc ) {
                    //comm.fill_agehisto_secondaryDengue( pooled_histo[ sim ][ id_next_checkpt ], compM );
                    comm.fill_agehisto_parity( pooled_histo[ sim ][ id_next_checkpt_inc ], compM );
                    ++id_next_checkpt_inc;
    
                    if ( id_next_checkpt_inc < n_checkpts_inc )
                        t_next_checkpt_inc = checkpts_inc[id_next_checkpt_inc];
                    else
                        t_next_checkpt_inc = Tmax + 1;
                    
                    comm.flush_events();
                    
                }
            }
            
            // demography update step
            if ( doDemography ) {
                if ( t < Tsample )
                    comm.do_demography_step(t, compM, false, false);
                else
                    comm.do_demography_step(t, compM, false, true);
            }
            
            // strain introduction update step
            comm.do_intro_step(t, intro_strains, compM);
        }
        compM.clear_incidence();
        comm.flush_events();
    }
    
    return { pooled_histo, seroprevalence };
}


std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> simulate_DENV_detailed_incidence(
                   const int& Nsim,
                   const int& seed,
                   const int& Tmax,
                   const int& Tsample,
                   const int& N,
                   const int& nstrains,
                   const double& beta,
                   const double& epsi,
                   const double& sigma,
                   const double& mean_degree,
                   const double& lD,
                   const double& gammaDD1,
                   const double& gammaDD2,
                   const double& intro_rate,
                   const double& seasonal_force,
                   const double& seasonal_peak,
                   const double& age_shape,
                   const double& age_scale,
                   const std::vector<int>& whichDs,
                   const std::vector<int>& checkpts ) {
    
    std::vector<double> trs = std::vector<double>( nstrains, beta );
    std::vector<double> epsis = std::vector<double>( nstrains, epsi );
    std::vector<double> sigs = std::vector<double>( nstrains, sigma );
    
    double lZ = 365.;
    
    double rhoDZ = 1.;
    double rhoZD = 1.;
    double rhoZD0 = 1.;
    double csiDZ = 1.;
    double csiZD = 1.;
    double gamma1D = gammaDD1;
    double gamma1Z = 1.;
    double gamma2D = gammaDD2;
    double gamma2Z = 1.;
    
    std::vector<int> endemic_strains = {};
    std::vector<int> all_strains = {};
    std::vector<int> intro_strains = {};
        
    // output flags
    
    std::mt19937_64 mymt = std::mt19937_64(1000);
    mymt.seed( seed );
    assign_mt(&mymt);
    
    CompartmentalModel compM = CompartmentalModel( nstrains, trs, sigs, epsis );
    compM.set_immunity_params( lD * 365., lZ, rhoDZ, rhoZD, rhoZD0, csiDZ, csiZD, gamma1D, gamma1Z, gamma2D, gamma2Z );
    
    for ( int i = 0; i < nstrains; i++ ) {
        if ( compM.is_endemic( i ) )
            endemic_strains.push_back( i );
        all_strains.push_back( i );
    }
    
    for ( int strain: whichDs )
        intro_strains.push_back( strain );
   
    // prepare output
    int n_checkpts_inc = static_cast<int>( checkpts.size() );
    assert( n_checkpts_inc > 0 );
    std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> pooled_histo;
    pooled_histo = std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>>( Nsim, std::vector<std::vector<std::vector<std::vector<int>>>>( n_checkpts_inc, std::vector<std::vector<std::vector<int>>>( 4,std::vector<std::vector<int>>( 101, std::vector<int>( 2, 0 ) ) ) ) );
   
    
    
    for ( int sim = 0; sim < Nsim; ++sim ) {
        
        int id_next_checkpt_inc  = 0;
        
        int t_next_checkpt_inc   = checkpts[ id_next_checkpt_inc ];

            
        Community comm = Community(N,
                                   intro_rate,
                                   mean_degree,
                                   seasonal_force,
                                   0.,
                                   seasonal_peak,
                                   0 );
        
        comm.initialize_population_ageonly( age_scale, age_shape );

        
        // initialize high levels of immunity to endemic strains in population
        for ( int strain: intro_strains )
            comm.initialize_immunity_strain( strain, compM );

        for ( int t = 0; t < Tmax; ++t ) {
                        
            if ( t == Tsample ) {
                compM.clear_incidence(); // clear incidence now (will save space on hard disk slightly)
            }
            
            if ( t < Tsample )
                comm.do_dynamics_step( t, compM );                        // epi dynamics
            else {
                comm.do_dynamics_step_output( t, compM );
                            
                
                // output: incidence
                if ( t == t_next_checkpt_inc ) {
                    //comm.fill_agehisto_secondaryDengue( pooled_histo[ sim ][ id_next_checkpt ], compM );
                    comm.fill_agehisto_parity_singleDvsAll( pooled_histo[ sim ][ id_next_checkpt_inc ], compM );
                    ++id_next_checkpt_inc;
    
                    if ( id_next_checkpt_inc < n_checkpts_inc )
                        t_next_checkpt_inc = checkpts[id_next_checkpt_inc];
                    else
                        t_next_checkpt_inc = Tmax + 1;
                    
                    comm.flush_events();
                    
                }
            }
            
            // demography update step
           
            if ( t < Tsample )
                comm.do_demography_step(t, compM, false, false);
            else
                comm.do_demography_step(t, compM, false, true);

            
            // strain introduction update step
            comm.do_intro_step(t, intro_strains, compM);
        }
        compM.clear_incidence();
        comm.flush_events();
    }
    
    return pooled_histo;
    
    
}


std::pair<std::vector<std::vector<std::vector<std::vector<int>>>>, std::vector<std::vector<std::vector<std::vector<int>>>>> simulate_DENV_sero(const int& Nsim,
                   const int& seed,
                   const int& Tmax,
                   const int& Tsample,
                   const int& N,
                   const int& nstrains,
                   const double& beta,
                   const double& epsi,
                   const double& sigma,
                   const double& mean_degree,
                   const double& lD,
                   const double& gammaDD1,
                   const double& gammaDD2,
                   const double& intro_rate,
                   const double& seasonal_force,
                   const double& seasonal_peak,
                   const double& age_shape,
                   const double& age_scale,
                   const double& sero_lam1,
                   const double& sero_lam2,
                   const double& sero_a0,
                   const std::vector<int>& checkpts_inc,
                   const std::vector<int>& checkpts_sero )  {
    
    std::vector<double> trs = std::vector<double>( nstrains, beta );
    std::vector<double> epsis = std::vector<double>( nstrains, epsi );
    std::vector<double> sigs = std::vector<double>( nstrains, sigma );
    
    double lZ = 365.;
    
    double rhoDZ = 1.;
    double rhoZD = 1.;
    double rhoZD0 = 1.;
    double csiDZ = 1.;
    double csiZD = 1.;
    double gamma1D = gammaDD1;
    double gamma1Z = 1.;
    double gamma2D = gammaDD2;
    double gamma2Z = 1.;
    
    std::vector<int> endemic_strains = {};
    std::vector<int> all_strains = {};
    std::vector<int> intro_strains = {};
        
    // output flags
    
    std::mt19937_64 mymt = std::mt19937_64(1000);
    mymt.seed( seed );
    assign_mt(&mymt);
    
    CompartmentalModel compM = CompartmentalModel( nstrains, trs, sigs, epsis );
    compM.set_immunity_params( lD * 365., lZ, rhoDZ, rhoZD, rhoZD0, csiDZ, csiZD, gamma1D, gamma1Z, gamma2D, gamma2Z );
    
    
    for ( int i = 0; i < nstrains; i++ ) {
        if ( compM.is_endemic( i ) )
            endemic_strains.push_back( i );
        all_strains.push_back( i );
    }
    
    for ( int strain: { 0, 1, 2, 3 } )
        intro_strains.push_back( strain );
   
    // prepare output
    int n_checkpts_inc = static_cast<int>( checkpts_inc.size() );
    assert( n_checkpts_inc > 0 );
    std::vector<std::vector<std::vector<std::vector<int>>>> pooled_histo;
    pooled_histo = std::vector<std::vector<std::vector<std::vector<int>>>>( Nsim, std::vector<std::vector<std::vector<int>>>( n_checkpts_inc, std::vector<std::vector<int>>( 4, std::vector<int>( 101, 0 ) ) ) );
    
    int n_checkpts_sero = static_cast<int>( checkpts_sero.size() );
    assert( n_checkpts_sero > 0 );
    std::vector<std::vector<std::vector<std::vector<int>>>> seroprevalence;
    seroprevalence = std::vector<std::vector<std::vector<std::vector<int>>>>( Nsim, std::vector<std::vector<std::vector<int>>>( n_checkpts_sero, std::vector<std::vector<int>>( 5, std::vector<int>( 101, 0 ) ) ) );
    
    
    for ( int sim = 0; sim < Nsim; ++sim ) {
        
       
        int id_next_checkpt_inc  = 0;
        int id_next_checkpt_sero = 0;
        
        int t_next_checkpt_inc   = checkpts_inc[ id_next_checkpt_inc   ];
        int t_next_checkpt_sero  = checkpts_sero[ id_next_checkpt_sero ];

            
        Community comm = Community(N,
                                   intro_rate,
                                   mean_degree,
                                   seasonal_force,
                                   0.,
                                   seasonal_peak,
                                   0 );
        
        // initialize age and immunity according to serology
        comm.initialize_population_serology( age_scale, age_shape, sero_lam1, sero_lam2, sero_a0, compM );


        for ( int t = 0; t < Tmax; ++t ) {
                        
            if ( t == Tsample ) {
                compM.clear_incidence(); // clear incidence now (will save space on hard disk slightly)
            }
            
            if ( t < Tsample )
                comm.do_dynamics_step( t, compM );                        // epi dynamics
            else {
                comm.do_dynamics_step_output( t, compM );
                
                // output: seroprevalence
                if ( t == t_next_checkpt_sero ) {
                    
                    comm.fill_seroprevalenceD( t, seroprevalence[ sim ][ id_next_checkpt_sero ], compM );
                    
                    ++id_next_checkpt_sero;
    
                    if ( id_next_checkpt_sero < n_checkpts_sero )
                        t_next_checkpt_sero = checkpts_sero[id_next_checkpt_sero];
                    else
                        t_next_checkpt_sero = Tmax + 1;
                    
                }
                
                
                // output: incidence
                if ( t == t_next_checkpt_inc ) {
                    //comm.fill_agehisto_secondaryDengue( pooled_histo[ sim ][ id_next_checkpt ], compM );
                    comm.fill_agehisto_parity( pooled_histo[ sim ][ id_next_checkpt_inc ], compM );
                    ++id_next_checkpt_inc;
    
                    if ( id_next_checkpt_inc < n_checkpts_inc )
                        t_next_checkpt_inc = checkpts_inc[id_next_checkpt_inc];
                    else
                        t_next_checkpt_inc = Tmax + 1;
                    
                    comm.flush_events();
                    
                }
            }
            
            // demography update step
          
            if ( t < Tsample )
                comm.do_demography_step(t, compM, false, false);
            else
                comm.do_demography_step(t, compM, false, true);

            
            // strain introduction update step
            comm.do_intro_step(t, intro_strains, compM);
        }
        compM.clear_incidence();
        comm.flush_events();
    }
    
    return { pooled_histo, seroprevalence };
}

// Generate a population with realistic age structure and a DENV immune landscape obtained
// from a sero-catalytic model. Simulates introduction of ZIKV and returns all infection events

std::vector<std::vector<EventInfectionO>> simulate_ZIKV_sero_detailed( const int& Nsim,
                      const int& seed,
                      const int& Tmax,
                      const int& TintroZ,
                      const int& N,
                      const int& nstrains,
                      const double& betaD,
                      const double& betaZ,
                      const double& epsiD,
                      const double& epsiZ,
                      const double& sigmaD,
                      const double& sigmaZ,
                      const double& mean_degree,
                      const double& lD,
                      const double& lZ,
                      const double& intro_rate,
                      const double& seasonal_forceD,
                      const double& seasonal_forceZ,
                      const double& seasonal_peakD,
                      const double& seasonal_peakZ,
                      const double& rhoDZ,
                      const double& rhoZD,
                      const double& rhoZD0,
                      const double& csiDZ,
                      const double& csiZD,
                      const double& gamma1D,
                      const double& gamma1Z,
                      const double& gamma2D,
                      const double& gamma2Z,
                      const double& age_shape,
                      const double& age_scale,
                      const double& sero_lam1,
                      const double& sero_lam2,
                      const double& sero_a0 )
{
    
    // set strain-specific parameters
    
    std::vector<double> trs   = {};
    std::vector<double> epsis = {};
    std::vector<double> sigs  = {};
    
    for ( int i = 0; i < nstrains; i++ ) {
        if ( i < nstrains - 1 ) {
            trs.push_back( betaD );
            epsis.push_back( epsiD );
            sigs.push_back( sigmaD );
        }
        else {
            trs.push_back( betaZ );
            epsis.push_back( epsiZ );
            sigs.push_back( sigmaZ );
        }
    }
    
    std::vector<int> endemic_strains = {};
    std::vector<int> all_strains = {};
    std::vector<int> intro_strains = {};
        
    for ( int i = 0; i < nstrains; i++ ) {
        if ( i < nstrains - 1 )
            endemic_strains.push_back( i );
        all_strains.push_back( i );
    }
    
    for ( int strain: { 0, 1, 2, 3 } )
        intro_strains.push_back( strain );
        
    // create comparmental model
        
    CompartmentalModel compM = CompartmentalModel( nstrains, trs, sigs, epsis );
    compM.set_immunity_params( lD * 365., lZ * 365., rhoDZ, rhoZD, rhoZD0, csiDZ, csiZD, gamma1D, gamma1Z, gamma2D, gamma2Z );
        
    // seed the rng
    
    std::mt19937_64 mymt = std::mt19937_64( 1000 );
    mymt.seed( seed );
    assign_mt(&mymt);
    
    // prepare vector with output
    std::vector<std::vector<EventInfectionO>> output = std::vector<std::vector<EventInfectionO>>( Nsim, std::vector<EventInfectionO>() );
    
    for ( int sim = 0; sim < Nsim; ++sim )
        output[sim].reserve( 20000 );
    
    
    // loop over simulations
    
    for ( int sim = 0; sim < Nsim; ++sim ) {
        
            
        Community comm = Community(N,
                                   intro_rate,
                                   mean_degree,
                                   seasonal_forceD,
                                   seasonal_forceZ,
                                   seasonal_peakD,
                                   seasonal_peakZ );
        
        // initialize age and immunity according to serology
        comm.initialize_population_serology( age_scale, age_shape, sero_lam1, sero_lam2, sero_a0, compM );


        for ( int t = 0; t < Tmax; ++t ) {
            
            // introduce ZIKV
            if ( t == TintroZ )
                comm.set_random_infection( compM.get_Z_id(), 10, t, compM );
            
           
            comm.do_dynamics_step_output( t, compM );
                
            comm.do_demography_step(t, compM, false, false);
            
            // strain introduction update step
            comm.do_intro_step( t, intro_strains, compM );
        }
        
        comm.swap_infection_events( &output[sim] );
        
        compM.clear_incidence();
        comm.flush_events();
    }
    
    return output;
    
}



// Simulate DENV. Remove a single strain for some time and then introduce it back

std::vector<std::vector<std::vector<std::vector<int>>>> simulate_DENV_removal(const int& Nsim,
                   const int& seed,
                   const int& Tmax,
                   const int& Tremoval,
                   const int& dT,
                   const int& N,
                   const int& nstrains,
                   const double& beta,
                   const double& epsi,
                   const double& sigma,
                   const double& mean_degree,
                   const double& lD,
                   const double& intro_rate,
                   const double& seasonal_force,
                   const double& seasonal_peak,
                   const double& age_shape,
                   const double& age_max_prob,
                   const std::vector<int>& whichDs,
                   const std::vector<int>& checkpts ) {
    
    std::vector<double> trs = std::vector<double>( nstrains, beta );
    std::vector<double> epsis = std::vector<double>( nstrains, epsi );
    std::vector<double> sigs = std::vector<double>( nstrains, sigma );

    double tr0 = trs[0];
    
    double lZ = 365.;

    double rhoDZ = 1.;
    double rhoZD = 1.;
    double rhoZD0 = 1.;
    double csiDZ = 1.;
    double csiZD = 1.;
    double gamma1D = 1.;
    double gamma1Z = 1.;
    double gamma2D = 1.;
    double gamma2Z = 1.;

    std::vector<int> endemic_strains = {};
    std::vector<int> all_strains = {};
    std::vector<int> intro_strains = {};
        
    // output flags

    std::mt19937_64 mymt = std::mt19937_64(1000);
    mymt.seed( seed );
    assign_mt(&mymt);

    CompartmentalModel compM = CompartmentalModel( nstrains, trs, sigs, epsis );
    compM.set_immunity_params( lD * 365., lZ, rhoDZ, rhoZD, rhoZD0, csiDZ, csiZD, gamma1D, gamma1Z, gamma2D, gamma2Z );

    for ( int i = 0; i < nstrains; i++ ) {
        if ( compM.is_endemic( i ) )
            endemic_strains.push_back( i );
        all_strains.push_back( i );
    }

    // prepare output
    
    int Tsample = Tremoval + dT;
    
    int n_checkpts = static_cast<int>( checkpts.size() );
    assert( n_checkpts > 0 );
    std::vector<std::vector<std::vector<std::vector<int>>>> pooled_histo;
    pooled_histo = std::vector<std::vector<std::vector<std::vector<int>>>>( Nsim, std::vector<std::vector<std::vector<int>>>( n_checkpts, std::vector<std::vector<int>>( 4, std::vector<int>( 101, 0 ) ) ) );


    for ( int sim = 0; sim < Nsim; ++sim ) {
        
        int t_next_checkpt = Tmax + 1;
        int id_next_checkpt = 0;
        
        if ( n_checkpts > 0 )
            t_next_checkpt = checkpts[id_next_checkpt];
        
        Community comm = Community(N,
                                   intro_rate,
                                   mean_degree,
                                   seasonal_force,
                                   0.,
                                   seasonal_peak,
                                   0 );
        
        comm.initialize_population_ageonly( age_max_prob, age_shape );

        
        // initialize high levels of immunity to endemic strains in population
        
        intro_strains.clear();
        for ( int strain: whichDs )
            intro_strains.push_back( strain );

        
        for ( int strain: whichDs )
            comm.initialize_immunity_strain( strain, compM );

        for ( int t = 0; t < Tmax; ++t ) {
            
            // facilitate disappearance of strain 0 by:
            // i)  Set its transmissibility to 0
            // ii) Prevent its re-introduction
            if ( t == Tremoval ) {
                std::remove( intro_strains.begin(), intro_strains.end(), 0 );
                compM.set_tr( 0, 0. );
            }
            
            // re-instate strain 0 back into the simulation
            if ( t == Tsample ) {
                compM.clear_incidence();
                intro_strains.push_back( 0 );
                compM.set_tr( 0, tr0 );
                comm.set_random_infection( 0, 5, t, compM );
            }
            
            if ( t < Tsample )
                comm.do_dynamics_step( t, compM );                        // epi dynamics
            else {
                comm.do_dynamics_step_output( t, compM );
                if ( t == t_next_checkpt ) {
                    //comm.fill_agehisto_secondaryDengue( pooled_histo[ sim ][ id_next_checkpt ], compM );
                    comm.fill_agehisto_parity( pooled_histo[ sim ][ id_next_checkpt ], compM );
                    ++id_next_checkpt;

                    if ( id_next_checkpt < n_checkpts )
                        t_next_checkpt = checkpts[id_next_checkpt];
                    else
                        t_next_checkpt = Tmax + 1;
                    
                    comm.flush_events();
                    
                }
            }
            
            // demography update step
            if ( t < Tsample )
                comm.do_demography_step(t, compM, false, false);
            else
                comm.do_demography_step(t, compM, false, true);
            
            // strain introduction update step
            comm.do_intro_step(t, intro_strains, compM);
        }
        compM.clear_incidence();
        comm.flush_events();
    }

    return pooled_histo;
}

// Simulate ZIKV and DENV

std::vector<std::vector<std::vector<std::vector<int>>>> simulate_ZIKV(   const int& Nsim,
                      const int& seed,
                      const int& Tmax,
                      const int& Tsample,
                      const int& TintroZ,
                      const int& N,
                      const int& nstrains,
                      const std::vector<int>& whichDs,
                      const std::vector<int>& checkPts,
                      const bool& printZonly,
                      const double& betaD,
                      const double& betaZ,
                      const double& epsiD,
                      const double& epsiZ,
                      const double& sigmaD,
                      const double& sigmaZ,
                      const double& mean_degree,
                      const double& lD,
                      const double& lZ,
                      const double& intro_rate,
                      const double& seasonal_forceD,
                      const double& seasonal_forceZ,
                      const double& seasonal_peakD,
                      const double& seasonal_peakZ,
                      const double& age_shape,
                      const double& age_max_prob,
                      const double& rhoDZ,
                      const double& rhoZD,
                      const double& rhoZD0,
                      const double& csiDZ,
                      const double& csiZD,
                      const double& gamma1D,
                      const double& gamma1Z,
                      const double& gamma2D,
                      const double& gamma2Z,
                      const int& nZIKV0 ) {
    
    // set strain-specific parameters
    
    std::vector<double> trs   = {};
    std::vector<double> epsis = {};
    std::vector<double> sigs  = {};
    
    for ( int i = 0; i < nstrains; i++ ) {
        if ( i < nstrains - 1 ) {
            trs.push_back( betaD );
            epsis.push_back( epsiD );
            sigs.push_back( sigmaD );
        }
        else {
            trs.push_back( betaZ );
            epsis.push_back( epsiZ );
            sigs.push_back( sigmaZ );
        }
    }
    
    std::vector<int> endemic_strains = {};
    std::vector<int> all_strains = {};
    std::vector<int> intro_strains = {};
        
    for ( int i = 0; i < nstrains; i++ ) {
        if ( i < nstrains - 1 )
            endemic_strains.push_back( i );
        all_strains.push_back( i );
    }
    
    for ( int strain: whichDs )
        intro_strains.push_back( strain );
        
    // create comparmental model
        
    CompartmentalModel compM = CompartmentalModel( nstrains, trs, sigs, epsis );
    compM.set_immunity_params( lD * 365., lZ * 365., rhoDZ, rhoZD, rhoZD0, csiDZ, csiZD, gamma1D, gamma1Z, gamma2D, gamma2Z );
        
    // seed the rng
    
    std::mt19937_64 mymt = std::mt19937_64( 1000 );
    mymt.seed( seed );
    assign_mt(&mymt);
    
    
    // prepare output vector
    // it is a ( Nsim * len(checkPts) * 3 * 101 ) vector
        
    assert( checkPts[0] > Tsample );
    int n_checkpts = static_cast<int>( checkPts.size() );
    
    int nhistos = 0;
    if ( printZonly )
        nhistos = 4;
    else
        nhistos = 12;
        
    std::vector<std::vector<std::vector<std::vector<int>>>> histo = std::vector<std::vector<std::vector<std::vector<int>>>>( Nsim ,std::vector<std::vector<std::vector<int>>>( checkPts.size(), std::vector<std::vector<int>>( nhistos, std::vector<int>( 101, 0 ) ) ) );
        
    // loop over simulations
    
    for ( int sim = 0; sim < Nsim; ++sim ) {
        
        int t_next_checkpt = Tmax + 1;
        int id_next_checkpt = 0;
        
        if ( n_checkpts > 0 )
            t_next_checkpt = checkPts[id_next_checkpt];
        
        compM.clear_incidence();
        Community comm = Community( N, intro_rate, mean_degree, seasonal_forceD, seasonal_forceZ, seasonal_peakD, seasonal_peakZ );
        comm.initialize_population_ageonly( age_max_prob, age_shape );

        
                
        // initialize population w/ high levels of immunity to DENV
        
        for ( int strain: whichDs )
            comm.initialize_immunity_strain( strain, compM );
        
        // loop over timesteps
        
        for ( int t = 0; t < Tmax; ++t ) {
            
            // check Zika introduction
            
            if ( t == TintroZ ) {
                comm.set_random_infection( compM.get_Z_id(), nZIKV0, t, compM);

                if ( TintroZ == Tsample ) {
                    compM.clear_incidence();
                }
                    
            }
            
            // run epi update
            
            if ( t < Tsample ) {
                comm.do_dynamics_step( t, compM ); // not recording infection events
            }
            else {
                comm.do_dynamics_step_output( t, compM ); // recording infection events
            }
                        
            // check output
            
            if ( t == t_next_checkpt ) {
                
                if ( printZonly )
                    comm.fill_agehisto_parity_Z( histo[sim][id_next_checkpt], compM );
                else
                    comm.fill_agehisto_parity_ZD( histo[sim][id_next_checkpt], compM );

                ++id_next_checkpt;

                if ( id_next_checkpt < n_checkpts )
                    t_next_checkpt = checkPts[id_next_checkpt];
                else
                    t_next_checkpt = Tmax + 1;
                
                comm.flush_events();
                
            }
            
            // run demography update
            
            comm.do_demography_step(t, compM, false, false);
            
            // strain introduction update step
            
            comm.do_intro_step(t, intro_strains, compM);

        }
        
        compM.clear_incidence();
        comm.flush_events();
    }
    return histo;
}


std::vector<std::vector<std::vector<std::vector<int>>>> simulate_ZIKV_sero(const int& Nsim,
                                                                      const int& seed,
                                                                      const int& Tmax,
                                                                      const int& Tsample,
                                                                      const int& TintroZ,
                                                                      const int& N,
                                                                      const int& nstrains,
                                                                      const std::vector<int>& whichDs,
                                                                      const std::vector<int>& checkPts,
                                                                      const bool& printZonly,
                                                                      const double& betaD,
                                                                      const double& betaZ,
                                                                      const double& epsiD,
                                                                      const double& epsiZ,
                                                                      const double& sigmaD,
                                                                      const double& sigmaZ,
                                                                      const double& mean_degree,
                                                                      const double& lD,
                                                                      const double& lZ,
                                                                      const double& intro_rate,
                                                                      const double& seasonal_forceD,
                                                                      const double& seasonal_forceZ,
                                                                      const double& seasonal_peakD,
                                                                      const double& seasonal_peakZ,
                                                                      const double& age_shape,
                                                                      const double& age_scale,
                                                                      const double& rhoDZ,
                                                                      const double& rhoZD,
                                                                      const double& rhoZD0,
                                                                      const double& csiDZ,
                                                                      const double& csiZD,
                                                                      const double& gamma1D,
                                                                      const double& gamma1Z,
                                                                      const double& gamma2D,
                                                                      const double& gamma2Z,
                                                                      const double& sero_lam1,
                                                                      const double& sero_lam2,
                                                                      const double& sero_a0,
                                                                      const int& nZIKV0 ) {
    
    
    // set strain-specific parameters
    
    std::vector<double> trs   = {};
    std::vector<double> epsis = {};
    std::vector<double> sigs  = {};
    
    for ( int i = 0; i < nstrains; i++ ) {
        if ( i < nstrains - 1 ) {
            trs.push_back( betaD );
            epsis.push_back( epsiD );
            sigs.push_back( sigmaD );
        }
        else {
            trs.push_back( betaZ );
            epsis.push_back( epsiZ );
            sigs.push_back( sigmaZ );
        }
    }
    
    std::vector<int> endemic_strains = {};
    std::vector<int> all_strains = {};
    std::vector<int> intro_strains = {};
        
    for ( int i = 0; i < nstrains; i++ ) {
        if ( i < nstrains - 1 )
            endemic_strains.push_back( i );
        all_strains.push_back( i );
    }
    
    for ( int strain: whichDs )
        intro_strains.push_back( strain );
        
    // create comparmental model
        
    CompartmentalModel compM = CompartmentalModel( nstrains, trs, sigs, epsis );
    compM.set_immunity_params( lD * 365., lZ * 365., rhoDZ, rhoZD, rhoZD0, csiDZ, csiZD, gamma1D, gamma1Z, gamma2D, gamma2Z );
        
    // seed the rng
    
    std::mt19937_64 mymt = std::mt19937_64( 1000 );
    mymt.seed( seed );
    assign_mt(&mymt);
    
    
    // prepare output vector
    // it is a ( Nsim * len(checkPts) * 3 * 101 ) vector
        
    assert( checkPts[0] > Tsample );
    int n_checkpts = static_cast<int>( checkPts.size() );
    
    int nhistos = 0;
    if ( printZonly )
        nhistos = 4;
    else
        nhistos = 12;
        
    std::vector<std::vector<std::vector<std::vector<int>>>> histo = std::vector<std::vector<std::vector<std::vector<int>>>>( Nsim, std::vector<std::vector<std::vector<int>>>( checkPts.size(), std::vector<std::vector<int>>( nhistos, std::vector<int>( 101, 0 ) ) ) );
        
    // loop over simulations
    
    for ( int sim = 0; sim < Nsim; ++sim ) {
        
        int t_next_checkpt = Tmax + 1;
        int id_next_checkpt = 0;
        
        if ( n_checkpts > 0 )
            t_next_checkpt = checkPts[id_next_checkpt];
        
        compM.clear_incidence();
        Community comm = Community( N, intro_rate, mean_degree, seasonal_forceD, seasonal_forceZ, seasonal_peakD, seasonal_peakZ );
        
        comm.initialize_population_serology( age_scale, age_shape, sero_lam1, sero_lam2, sero_a0, compM );
        
      
        
        // loop over timesteps
        
        for ( int t = 0; t < Tmax; ++t ) {
            
            // check Zika introduction
            
            if ( t == TintroZ ) {
                comm.set_random_infection( compM.get_Z_id(), nZIKV0, t, compM);

                if ( TintroZ == Tsample ) {
                    compM.clear_incidence();
                }
                    
            }
            
            // run epi update
            
            if ( t < Tsample ) {
                comm.do_dynamics_step( t, compM ); // not recording infection events
            }
            else {
                comm.do_dynamics_step_output( t, compM ); // recording infection events
            }
                        
            // check output
            
            if ( t == t_next_checkpt ) {
                
                if ( printZonly )
                    comm.fill_agehisto_parity_Z( histo[sim][id_next_checkpt], compM );
                else
                    comm.fill_agehisto_parity_ZD( histo[sim][id_next_checkpt], compM );

                ++id_next_checkpt;

                if ( id_next_checkpt < n_checkpts )
                    t_next_checkpt = checkPts[id_next_checkpt];
                else
                    t_next_checkpt = Tmax + 1;
                
                comm.flush_events();
                
            }
            
            // run demography update
            
            comm.do_demography_step(t, compM, false, false);
            
            // strain introduction update step
            
            comm.do_intro_step(t, intro_strains, compM);

        }
        
        compM.clear_incidence();
        comm.flush_events();
    }
    return histo;
    
}


/*
 Simulate the dynamics of ZIKV
 Uses realistic seroprevalence profile
 Measures serology
 */

/*
std::pair<std::vector<std::vector<std::vector<std::vector<int>>>>, std::vector<std::vector<std::vector<std::vector<int>>>>> simulate_ZIKV_sero_getSerology(const int& Nsim,
                                                                      const int& seed,
                                                                      const int& Tmax,
                                                                      const int& Tsample,
                                                                      const int& TintroZ,
                                                                      const int& N,
                                                                      const int& nstrains,
                                                                      const std::vector<int>& whichDs,
                                                                      const std::vector<int>& checkpts_inc,
                                                                      const std::vector<int>& checkpts_sero,
                                                                      const bool& printZonly,
                                                                      const double& betaD,
                                                                      const double& betaZ,
                                                                      const double& epsiD,
                                                                      const double& epsiZ,
                                                                      const double& sigmaD,
                                                                      const double& sigmaZ,
                                                                      const double& mean_degree,
                                                                      const double& lD,
                                                                      const double& lZ,
                                                                      const double& intro_rate,
                                                                      const double& seasonal_forceD,
                                                                      const double& seasonal_forceZ,
                                                                      const double& seasonal_peakD,
                                                                      const double& seasonal_peakZ,
                                                                      const double& age_shape,
                                                                      const double& age_max_prob,
                                                                      const double& rhoDZ,
                                                                      const double& rhoZD,
                                                                      const double& rhoZD0,
                                                                      const double& csiDZ,
                                                                      const double& csiZD,
                                                                      const double& gamma1D,
                                                                      const double& gamma1Z,
                                                                      const double& gamma2D,
                                                                      const double& gamma2Z,
                                                                      const double& sero_lam1,
                                                                      const double& sero_lam2,
                                                                      const double& sero_a0,
                                                                      const int& nZIKV0 = 10 )  {
    
    
    // set strain-specific parameters
    
    std::vector<double> trs   = {};
    std::vector<double> epsis = {};
    std::vector<double> sigs  = {};
    
    for ( int i = 0; i < nstrains; i++ ) {
        if ( i < nstrains - 1 ) {
            trs.push_back( betaD );
            epsis.push_back( epsiD );
            sigs.push_back( sigmaD );
        }
        else {
            trs.push_back( betaZ );
            epsis.push_back( epsiZ );
            sigs.push_back( sigmaZ );
        }
    }
    
    std::vector<int> endemic_strains = {};
    std::vector<int> all_strains = {};
    std::vector<int> intro_strains = {};
        
    for ( int i = 0; i < nstrains; i++ ) {
        if ( i < nstrains - 1 )
            endemic_strains.push_back( i );
        all_strains.push_back( i );
    }
    
    for ( int strain: whichDs )
        intro_strains.push_back( strain );
        
    // create comparmental model
        
    CompartmentalModel compM = CompartmentalModel( nstrains, trs, sigs, epsis );
    compM.set_immunity_params( lD * 365., lZ * 365., rhoDZ, rhoZD, rhoZD0, csiDZ, csiZD, gamma1D, gamma1Z, gamma2D, gamma2Z );
        
    // seed the rng
    
    std::mt19937_64 mymt = std::mt19937_64( 1000 );
    mymt.seed( seed );
    assign_mt(&mymt);
    
    
    // prepare output

    int nhistos = 0;
    if ( printZonly )
        nhistos = 4;
    else
        nhistos = 12;
    
    int n_checkpts_inc = static_cast<int>( checkpts_inc.size() );
    assert( n_checkpts_inc > 0 );
    std::vector<std::vector<std::vector<std::vector<int>>>> pooled_histo;
    pooled_histo = std::vector<std::vector<std::vector<std::vector<int>>>>( Nsim, std::vector<std::vector<std::vector<int>>>( n_checkpts_inc, std::vector<std::vector<int>>( nhistos, std::vector<int>( 101, 0 ) ) ) );
    
    int n_checkpts_sero = static_cast<int>( checkpts_sero.size() );
    assert( n_checkpts_sero > 0 );
    std::vector<std::vector<std::vector<std::vector<int>>>> seroprevalence;
    seroprevalence = std::vector<std::vector<std::vector<std::vector<int>>>>( Nsim, std::vector<std::vector<std::vector<int>>>( n_checkpts_sero, std::vector<std::vector<int>>( 5, std::vector<int>( 101, 0 ) ) ) );
    
        
    // loop over simulations
    
    for ( int sim = 0; sim < Nsim; ++sim ) {
       
        
        int id_next_checkpt_inc  = 0;
        int id_next_checkpt_sero = 0;
        
        int t_next_checkpt_inc   = checkpts_inc[ id_next_checkpt_inc   ];
        int t_next_checkpt_sero  = checkpts_sero[ id_next_checkpt_sero ];
        
        
        
        compM.clear_incidence();
        Community comm = Community( N, intro_rate, mean_degree, seasonal_forceD, seasonal_forceZ, seasonal_peakD, seasonal_peakZ );
        
        comm.initialize_population_serology( age_scale, age_shape, sero_lam1, sero_lam2, sero_a0, compM );
        
      
        for ( int t = 0; t < Tmax; ++t ) {
            
            // check Zika introduction
            
            if ( t == TintroZ ) {
                comm.set_random_infection( compM.get_Z_id(), nZIKV0, t, compM);

                if ( TintroZ == Tsample ) {
                    compM.clear_incidence();
                }
                    
            }
            
            // run epi update
            
            if ( t < Tsample ) {
                comm.do_dynamics_step( t, compM ); // not recording infection events
            }
            else {
                comm.do_dynamics_step_output( t, compM ); // recording infection events
            }
                        
            // check output
            
            // !!! remove this section
            if ( t == t_next_checkpt ) {
                
                if ( printZonly )
                    comm.fill_agehisto_parity_Z( histo[sim][id_next_checkpt], compM );
                else
                    comm.fill_agehisto_parity_ZD( histo[sim][id_next_checkpt], compM );

                ++id_next_checkpt;

                if ( id_next_checkpt < n_checkpts )
                    t_next_checkpt = checkPts[id_next_checkpt];
                else
                    t_next_checkpt = Tmax + 1;
                
                comm.flush_events();
                
            }
            
            // output: seroprevalence
            if ( t == t_next_checkpt_sero ) {
                
                comm.fill_seroprevalenceD( t, seroprevalence[ sim ][ id_next_checkpt_sero ], compM );
                
                ++id_next_checkpt_sero;

                if ( id_next_checkpt_sero < n_checkpts_sero )
                    t_next_checkpt_sero = checkpts_sero[id_next_checkpt_sero];
                else
                    t_next_checkpt_sero = Tmax + 1;
                
            }
            
            
            // output: incidence
            if ( t == t_next_checkpt_inc ) {
                //comm.fill_agehisto_secondaryDengue( pooled_histo[ sim ][ id_next_checkpt ], compM );
                comm.fill_agehisto_parity( pooled_histo[ sim ][ id_next_checkpt_inc ], compM );
                ++id_next_checkpt_inc;

                if ( id_next_checkpt_inc < n_checkpts_inc )
                    t_next_checkpt_inc = checkpts_inc[id_next_checkpt_inc];
                else
                    t_next_checkpt_inc = Tmax + 1;
                
                comm.flush_events();
                
            }
            
            
            // run demography update
            
            comm.do_demography_step(t, compM, false, false);
        
            
            // strain introduction update step
            
            comm.do_intro_step(t, intro_strains, compM);

        }
        
        compM.clear_incidence();
        comm.flush_events();
    }
    return { pooled_histo, seroprevalence } ;
    
}
*/






std::vector<std::vector<std::vector<std::vector<int>>>> simulate_ZIKV_parityDependentCrossReactivity(const int& Nsim,
                                                  const int& seed,
                                                  const int& Tmax,
                                                  const int& Tsample,
                                                  const int& TintroZ,
                                                  const int& N,
                                                  const int& nstrains,
                                                  const std::vector<int>& whichDs,
                                                  const std::vector<int>& checkPts,
                                                  const bool& printZonly,
                                                  const double& betaD,
                                                  const double& betaZ,
                                                  const double& epsiD,
                                                  const double& epsiZ,
                                                  const double& sigmaD,
                                                  const double& sigmaZ,
                                                  const double& mean_degree,
                                                  const double& lD,
                                                  const double& lZ,
                                                  const double& intro_rate,
                                                  const double& seasonal_forceD,
                                                  const double& seasonal_forceZ,
                                                  const double& seasonal_peakD,
                                                  const double& seasonal_peakZ,
                                                  const double& age_shape,
                                                  const double& age_max_prob,
                                                  const double& rhoZD0,
                                                  const double& rhoZD1,
                                                  const double& rhoZD2,
                                                  const double& chiZD0,
                                                  const double& chiZD1,
                                                  const double& chiZD2 ) {
    
    
    // set strain-specific parameters
    
    std::vector<double> trs   = {};
    std::vector<double> epsis = {};
    std::vector<double> sigs  = {};
    
    for ( int i = 0; i < nstrains; i++ ) {
        if ( i < nstrains - 1 ) {
            trs.push_back( betaD );
            epsis.push_back( epsiD );
            sigs.push_back( sigmaD );
        }
        else {
            trs.push_back( betaZ );
            epsis.push_back( epsiZ );
            sigs.push_back( sigmaZ );
        }
    }
    
    std::vector<int> endemic_strains = {};
    std::vector<int> all_strains = {};
    std::vector<int> intro_strains = {};
        
    for ( int i = 0; i < nstrains; i++ ) {
        if ( i < nstrains - 1 )
            endemic_strains.push_back( i );
        all_strains.push_back( i );
    }
    
    for ( int strain: whichDs )
        intro_strains.push_back( strain );
        
    // create comparmental model
        
    CompartmentalModel compM = CompartmentalModel( nstrains, trs, sigs, epsis );
    compM.set_immunity_params_v2( lD * 365., lZ * 365., rhoZD0, rhoZD1, rhoZD2, chiZD0, chiZD1, chiZD2 );
        
    // seed the rng
    
    std::mt19937_64 mymt = std::mt19937_64( 1000 );
    mymt.seed( seed );
    assign_mt(&mymt);
    
    
    // prepare output vector (incidence)
    // it is a ( Nsim * len(checkPts) * 3 * 101 ) vector
        
    assert( checkPts[0] > Tsample );
    int n_checkpts = static_cast<int>( checkPts.size() );
    
    int nhistos = 0;
    if ( printZonly )
        nhistos = 3;
    else
        nhistos = 9;
        
    std::vector<std::vector<std::vector<std::vector<int>>>> histo = std::vector<std::vector<std::vector<std::vector<int>>>>( Nsim ,std::vector<std::vector<std::vector<int>>>( checkPts.size(), std::vector<std::vector<int>>( nhistos, std::vector<int>( 101, 0 ) ) ) );
    
    
    
    // loop over simulations
    
    for ( int sim = 0; sim < Nsim; ++sim ) {
        
        int t_next_checkpt = Tmax + 1;
        int id_next_checkpt = 0;
        
        if ( n_checkpts > 0 )
            t_next_checkpt = checkPts[id_next_checkpt];
        
        compM.clear_incidence();
        Community comm = Community( N, intro_rate, mean_degree, seasonal_forceD, seasonal_forceZ, seasonal_peakD, seasonal_peakZ );
        comm.initialize_population_ageonly( age_max_prob, age_shape );

        
                
        // initialize population w/ high levels of immunity to DENV
        
        for ( int strain: whichDs )
            comm.initialize_immunity_strain( strain, compM );
        
        // loop over timesteps
        
        for ( int t = 0; t < Tmax; ++t ) {
            
            // check Zika introduction
            
            if ( t == TintroZ ) {
                
                comm.set_random_infection( compM.get_Z_id(), 25, t, compM);

                if ( TintroZ == Tsample ) {
                    
                    compM.clear_incidence();
                
                }
                    
            }
            
            // run epi update
            
            if ( t < Tsample ) {
                
                comm.do_dynamics_step( t, compM ); // not recording infection events
            
            }
            else {
            
                comm.do_dynamics_step_output( t, compM ); // recording infection events
            
            }
                        
            // check output
            
            if ( t == t_next_checkpt ) {
                
                if ( printZonly )
                    comm.fill_agehisto_parity_Z( histo[sim][id_next_checkpt], compM );
                else
                    comm.fill_agehisto_parity_ZD( histo[sim][id_next_checkpt], compM );

                ++id_next_checkpt;

                if ( id_next_checkpt < n_checkpts )
                    t_next_checkpt = checkPts[id_next_checkpt];
                else
                    t_next_checkpt = Tmax + 1;
                
                comm.flush_events();
                
            }
            
            // run demography update
            
            comm.do_demography_step(t, compM, false, false);
            
            // strain introduction update step
            
            comm.do_intro_step(t, intro_strains, compM);

        }
        
        compM.clear_incidence();
        comm.flush_events();
    }
    return histo;
}
