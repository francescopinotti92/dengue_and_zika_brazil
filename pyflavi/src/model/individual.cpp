//
//  individual.cpp
//  CrossFlavivirus
//
//  Created by MacBook Pro on 05/01/2021.
//

#include "individual.hpp"

mt19937_64* mt = nullptr;
uniform_int_distribution<int> randint_population = uniform_int_distribution<int>(0, 1);
uniform_real_distribution<double> uniGen = uniform_real_distribution<double>(0., 1.);
uniform_int_distribution<int> randint_dayofyear = uniform_int_distribution<int>(0, 364);
exponential_distribution<double> expoGen = exponential_distribution<double>(1.);

void assign_mt(mt19937_64* mt_) {
    mt = mt_;
};

AgeSamplerClass::AgeSamplerClass(int age_max_) {
    age_max = age_max_;
    psurv.resize(0);
    psurv.resize(age_max, 0.);

    for(int a = 0; a < age_max; ++a) {
        double tmp = 1.;
        for(int a1 = 0; a1 <= a; ++a1)
            tmp *= 1. - ageDependentMortality(a1 * 1.0);
        psurv[a] = tmp;
    }
};

double ageDependentMortality( const double& age ){
    //return( 1. - exp( -pow(age/MAXAGE_DEATHPROB,SHAPE_DEATHPROB) ) );
    return( pow( (age + 1.)/MAXAGE_DEATHPROB,SHAPE_DEATHPROB) - pow(age/MAXAGE_DEATHPROB,SHAPE_DEATHPROB) );
};

int sampleAge() {
    static AgeSamplerClass age_distr(100);
    
    int age = uniGen(*mt) * age_distr.age_max;
    while(true) {
        if (uniGen(*mt) < age_distr.psurv[age])
            break;
        else
            age = uniGen(*mt) * age_distr.age_max;
    }
    return (int)(age * 365); // convert in days
};



//====== Compartmental model class ======//

CompartmentalModel::CompartmentalModel(int nstrains_, vector<double> tr, vector<double> sig, vector<double> ltcy): nstrains(nstrains_), tr_i(tr), sig_i(sig), ltc_prob_i(ltcy) {
    
    // compute recovery probability
    rec_prob_i.resize(0);
    for (int i = 0; i < nstrains; ++i) {
        //rec_prob_i.push_back(1. - exp( -sig_i[i] ) ); // transform parameters
        rec_prob_i.push_back( sig_i[i] ); // no transformation
    }
    
    strains_carried.resize( (int)pow( 2, nstrains ), {}); // do not resize this!!
    strains_carried_bool.resize( (int)pow(2 , nstrains ), vector<bool>( nstrains, false ) );
    
    // Enumerate combinations of strains for each infectious state
    for ( int i = 0; i < (int)pow(2, nstrains); ++i ) {
        std::bitset<100> inf_state(i);
        for ( uint j = 0; j < nstrains; ++j ) {
            if (inf_state[j] == 1) {
                strains_carried[i].push_back( j ); // fill list of strains
                strains_carried_bool[i][j] = true;
            }
        }
    }
    
    cum_incidence.resize( nstrains, 0 );
    
    for (int i = 0; i < nstrains; ++i) {
        waningGens.push_back( geometric_distribution<int>(1.) );
    }
    
    for (int i = 0; i < nstrains; ++i) {
        latencyGens.push_back( geometric_distribution<int>( ltc_prob_i[i]) );
    }
    
    // default parameters
    
    gamma_noZ_D = vector<double>( nstrains - 1, 1. );
    gamma_ltZ_D = vector<double>( nstrains - 1, 1. );
    gamma_stZ_D = vector<double>( nstrains - 1, 1. );
    gamma_ltD_Z = vector<double>( nstrains, 1. );
    gamma_stD_Z = vector<double>( nstrains, 1. );
    
    csiDZ = 1.;
    csiZD = 1.;
    
}

void CompartmentalModel::set_immunity_params(double lD, double lZ, double rhoDZ, double rhoZD) {
   
    for (int i = 0; i < nstrains - 1; ++i) {
        waningGens[i].param( geometric_distribution<int>::param_type( 1 - exp(-1./lD) ) );
    }
    waningGens[nstrains - 1].param( geometric_distribution<int>::param_type( 1 - exp(-1./lZ) ) );

    
    gamma_stD_Z = vector<double>( nstrains - 1, rhoDZ );
    gamma_stZ_D = vector<double>( nstrains - 1, rhoZD );

    

}

// Only short term effects of ZIKV cross-reactivity to DENV
void CompartmentalModel::set_immunity_params(double lD, double lZ, double rhoDZ, double rhoZD, double rhoZD0, double csiDZ, double csiZD, double gamma1D, double gamma1Z, double gamma2D, double gamma2Z) {
    
    
    for (int i = 0; i < nstrains - 1; ++i) {
        waningGens[i].param( geometric_distribution<int>::param_type( 1 - exp(-1./lD) ) );
    }
    waningGens[nstrains - 1].param( geometric_distribution<int>::param_type( 1 - exp(-1./lZ) ) );

    
    this->csiDZ = csiDZ;
    this->csiZD = csiZD;
    
    gamma_noZ_D = { 1., gamma1D, gamma2D, gamma2D };
    gamma_ltD_Z = { 1., gamma1Z, gamma2Z, gamma2Z, gamma2Z };
    gamma_stD_Z = vector<double>( nstrains, rhoDZ );
    gamma_stZ_D = { rhoZD0, rhoZD, rhoZD, rhoZD};
    
    gamma_ltZ_D = vector<double>( nstrains - 1, 1. );
    
}

// short and long terms of ZIKV cross-reactivity to DENV, dependening on parity of DENV
void CompartmentalModel::set_immunity_params_v2( double lD, double lZ, double rhoZD0, double rhoZD1, double rhoZD2, double chiZD0, double chiZD1, double chiZD2 ) {
    
    // set distributions for waning of immunity
    for (int i = 0; i < nstrains - 1; ++i) {
        waningGens[i].param( geometric_distribution<int>::param_type( 1 - exp(-1./lD) ) );
    }
    waningGens[nstrains - 1].param( geometric_distribution<int>::param_type( 1 - exp(-1./lZ) ) );

    gamma_noZ_D = vector<double>( nstrains - 1, 1. );
    gamma_stZ_D = { rhoZD0, rhoZD1, rhoZD2, rhoZD2 };
    gamma_ltZ_D = { chiZD0, chiZD1, chiZD2, chiZD2 };
    
    gamma_ltD_Z = vector<double>( nstrains, 1. );
    gamma_stD_Z = vector<double>( nstrains, 1. );


    
}

double CompartmentalModel::get_tr( Individual& source, const int &strain ) {
    
    if ( is_endemic( strain ) ) {
        if ( source.get_ctrZ() > 0 )
            return csiZD * tr_i[strain];
        else
            return tr_i[strain];
    }
    else {
        if ( source.get_ctrD() > 0 )
            return csiDZ * tr_i[strain];
        else
            return tr_i[strain];
    }
    
}


/* Compute susceptibility of 'target' to 'strain' */
/*
double CompartmentalModel::get_immunity(Individual &target, const int& strain) {
    
    // i) target has already been exposed to strain => susc = 0
    
    if ( is_immune( strain, target ) )
        return 0.;
    
    // ii) target is infectious or latent => susc = 0
    
    if ( target.is_infectious() or target.is_latent() )
        return 0.; // forbid co-infections for simplicity

    // iii) check immune history
    
    if (  is_endemic( strain ) ) { // D infection
        
        // iii-a) short-term cross protection from D => susc = 0
        
        if ( target.get_ctrD() > 0 ) {
            return 0.; // cross_protected by D against D
        }
        else { // iii-b)
            double res = 1.;
            if ( target.get_ctrZ() > 0 ) { // check cross_imm from Zika
                if ( target.get_nD() == 0 )
                    res *= rhoZD0;
                else
                    res *= rhoZD;
            }
            
            if ( target.get_nD() == 1 ) // check 1D status
                res *= gamma1D;
            else if ( target.get_nD() >= 2 ) // check 2+D status
                res *= gamma2D;
            
            return res;
        }
    }
    else {  // Z infection
        if ( target.get_ctrZ() > 0 ) {
            return 0.; // cross_protected by Z against Z; ???
        }
        else {
            double res = 1.;
            if ( target.get_ctrD() > 0 ) // short-term cross-imm from D
                res *= rhoDZ;
            
            if ( target.get_nD() == 1 )      // check 1D status
                res *= gamma1Z;
            else if ( target.get_nD() >= 2 ) // check 2+D status
                res *= gamma2Z;
            return res;
        }
    }
    
    //return get_tr(strain) * cross_imm_susc[target.get_imm_state()][strain];
}
 */

/* Compute susceptibility of 'target' to 'strain' */
// modified version
double CompartmentalModel::get_immunity( Individual &target, const int& strain ) {
    
    // i) target has already been exposed to strain => susc = 0
    
    if ( is_immune( strain, target ) )
        return 0.;
    
    // ii) target is infectious or latent => susc = 0
    
    if ( target.is_infectious() or target.is_latent() )
        return 0.; // forbid co-infections for simplicity

    // iii) check immune history
    
    if (  is_endemic( strain ) ) { // D infection
        
        // iii-a) short-term cross protection from D => susc = 0
        
        if ( target.get_ctrD() > 0 ) {
            return 0.; // cross_protected by D against D
        }
        else { // iii-b)
            int nD = target.get_nD();
            int nZ = target.get_nZ();
            
            if ( target.get_ctrZ() == 0 ) {
                if ( nZ == 0 )
                    return gamma_noZ_D[nD];
                else
                    return gamma_ltZ_D[nD];
            }
            else
                return gamma_stZ_D[nD];
        }
    }
    else {  // Z infection
        if ( target.get_ctrZ() > 0 ) {
            return 0.; // cross_protected by Z against Z; ???
        }
        else {
            int nD = target.get_nD();
            if ( target.get_ctrD() > 0 ) // short-term cross-imm from D
                return gamma_stD_Z[nD];
            else {
                return gamma_ltD_Z[nD];
            }
        
        }
    }
    //return get_tr(strain) * cross_imm_susc[target.get_imm_state()][strain];
}


inline vector<int>& CompartmentalModel::get_carried_strains(Individual& indiv) { return strains_carried[indiv.get_inf_state()]; }

inline vector<int>& CompartmentalModel::get_imm_strains(Individual& indiv) { return strains_carried[indiv.get_imm_state()]; }

inline vector<int>& CompartmentalModel::get_carried_strains(const int& inf_state) { return strains_carried[inf_state]; }

inline bool CompartmentalModel::is_immune( const int& strain, const Individual& indiv ) {
    return strains_carried_bool[ indiv.get_imm_state() ][ strain ];
}


// shuffle order of strains in each state
void CompartmentalModel::shuffle_strain_set() {
    for (int i = 1, ns = (int)pow(2, nstrains); i < ns; ++i) {
        shuffle(strains_carried[i].begin(), strains_carried[i].end(), *mt);
    }
}

void CompartmentalModel::print_cumulative_incidence(int t, ostream& outfile) {
    outfile << t << "," << stringify_vector(cum_incidence) << endl;
}


//====== Individual class ======//

Individual::Individual(int id_, int index_, int dayborn_): id(id_), index(index_), dayborn(dayborn_) {
    inf_state = 0;
    imm_state = 0;
    imm_ctrD = 0;
    imm_ctrZ = 0;
    nD = 0;
    nZ = 0;
    tLastD  = numeric_limits<int>::min();
    tFirstD = numeric_limits<int>::min();
    tZ = numeric_limits<int>::min();
    hazard = expoGen(*mt);
    infection_evs.clear();
    infection_evs.reserve(5);
    was_updated = false;
    isL = false;
}

void Individual::death(int t, int newID) {
    dayborn = t;
    id = newID;
    inf_state = 0;
    imm_state = 0;
    imm_ctrD = 0;
    imm_ctrZ = 0;
    nD = 0;
    nZ = 0;
    tLastD  = numeric_limits<int>::min();
    tFirstD = numeric_limits<int>::min();
    tZ = numeric_limits<int>::min();
    hazard = expoGen(*mt);
    infection_evs.resize(0);
    was_updated = false;
    isL = false;
}

// returns true if hazard hits 0
bool Individual::update_hazard( double hz ) {
    
    if ( was_updated)
        return false; // only one update per timestep
    
    hazard -= hz;
    if (hazard <= 0.) {
        hazard = expoGen(*mt);
        was_updated = true;
        return true;
    }
    else
        return false;
    
    // try the standard method
    /*
    double pp = 1. - exp(-hz);
    if (uniGen(*mt) < pp) {
        was_updated = true;
        return true;
    }
    else
        return false;
    */
}

// returns true if was susceptible and now is infectious
bool Individual::apply_infection(const int& t, const int& strain, const bool& isD) {
    
    bool was_susceptible = !is_infectious();
    inf_state += (1 << strain);
    imm_state += (1 << strain);
    was_updated = false;
    add_infection_ev( strain, t );
    
    if ( isD )  {
        ++nD;
        tLastD = t;
        if ( tFirstD == numeric_limits<int>::min() ) {
            tFirstD = t;
        }
    }
    else {
        ++nZ;
        tZ = t;
    }
    
    set_latent_state(false);
    
    return was_susceptible;
}



// returns true if was infectious and now it isn't
bool Individual::apply_recovery( const int& strain ) {
    
    inf_state -= ( 1 << strain );
        
    if ( is_infectious() )
        return false;
    else
        return true;
}

void Individual::add_immunity( const int& strain ) {
    imm_state += ( 1 << strain );
}

void Individual::add_immunity( const int& strain, const int& t, const bool& isD ) {
    
    add_immunity( strain );
    
    if ( isD )  {
        ++nD;
        if ( t > tLastD )
            tLastD = t;
        if ( ( tFirstD == numeric_limits<int>::min() ) or ( t < tFirstD ) )
            tFirstD = t;
    }
    else {
        ++nZ;
        tZ = t;
    }
    
    
}


// check if individual is already infected with 'strain'
bool Individual::is_infectious_with(const int& strain, CompartmentalModel &compM) {
    vector<int>& strains = compM.get_carried_strains(*this);
    if ( find( strains.begin(), strains.end(), strain ) == strains.end() )
        return false;
    else
        return true;
}

void Individual::update_crossImm(const int& strain, const int& amount, CompartmentalModel &compM) {
    if ( compM.is_endemic(strain) )
        imm_ctrD += amount;
    else
        imm_ctrZ += amount;
}


// check if individual is already exposed to 'strain'
bool Individual::is_immune_to(const int& strain, CompartmentalModel& compM) {
    
    vector<int>& strains = compM.get_imm_strains(*this); // immune and infected states are the same in terms of structure
    if ( find(strains.begin(), strains.end(), strain) == strains.end() )
        return false;
    else
        return true;
}

// count past exposures to Dengue
int Individual::is_exposed_toD(CompartmentalModel& compM) {
    vector<int>& strains = compM.get_imm_strains(*this); // immune and infected states are the same in terms of structure
    int ct = 0;
    for ( int& strain: strains ) {
        if ( compM.is_endemic( strain ) )
            ct++;
    }
    return ct;
}

int Individual::is_exposed_toZ(CompartmentalModel &compM) {
    /*
    vector<int>& strains = compM.get_imm_strains( *this ); // immune and infected states are the same in terms of structure
    int ct = 0;
    for ( int& strain: strains ) {
        if ( compM.is_invasive( strain ) )
            ct++;
    }
    return ct;
    */
    return nZ; // this is a simple improvement
}


//====== Community class ======//

Community::Community( int N_, double intro_rate_, double mean_degree_, double seasonal_force_D_, double seasonal_force_Z_, double seasonal_max_t_D_, double seasonal_max_t_Z_ ): N(N_), inf_evs(N), rec_evs(N), intro_rate(intro_rate_), mean_degree(mean_degree_), seasonal_force_D(seasonal_force_D_), seasonal_force_Z(seasonal_force_Z_), seasonal_max_t_D( seasonal_max_t_D_ ), seasonal_max_t_Z( seasonal_max_t_Z_ ), age_manager( 3., 62. , 100 ) {
    
    nodes.resize(0);
    nodes.reserve(N);
    infected_nodes.clear();
    
    maxID = 0;
    
    //auto cmp = [](const DeathEvent& left, const DeathEvent& right) { return (left.t) > (right.t); };
    //death_events;
    
    /*
    for (unsigned int i = 0; i < N; ++i) {
        //nodes.push_back( Individual(maxID, i, -sampleAge() ) ); // it is negative because we begin at t = 0.
        
        pair<int, int> age_data = age_manager.sample_joint_age_residual_eq(*mt);
        nodes.push_back( Individual(maxID, i, - age_data.first ) );
        death_events.push( DeathEvent( i, age_data.second ) );
        maxID++;
    }*/

    randint_population.param( uniform_int_distribution<int>::param_type(0, N - 1) );
    poi_intros.param( poisson_distribution<int>::param_type(intro_rate) );
    poi_contacts.param( poisson_distribution<int>::param_type(mean_degree) );
    
    inf_events_output.reserve(3 * N);
    
    // allocate waning events
    vector<WaneEvent> container;
    container.reserve(N * 3);
    waning_events = priority_queue<WaneEvent, vector<WaneEvent>, CompareWaneEvents>( CompareWaneEvents(), std::move(container));
    
    // allocate latency events
    vector<InfEvent> containerL;
    containerL.reserve(N * 3);
    latency_events = priority_queue<InfEvent, vector<InfEvent>, CompareInfEvents>( CompareInfEvents(), std::move(containerL));
}

void Community::initialize_population_ageonly( double age_scale, double age_shape ) {
    
    // update age_manager if necessary
    age_manager = WeibullAgeSampler( age_shape, age_scale, 100 );
    
    nodes.clear();
    maxID = 0;

    while( !death_events.empty() )
        death_events.pop();
    
    for ( int i = 0; i < N; ++i ) {
        pair<int, int> age_data = age_manager.sample_joint_age_residual_eq(*mt);
        nodes.push_back( Individual( maxID, i, - age_data.first ) );
        death_events.push( DeathEvent( i, age_data.second ) );
        maxID++;
    }
    
}

void Community::initialize_population_serology( double age_scale, double age_shape, double sero_lam1, double sero_lam2, double sero_a0, CompartmentalModel& compM ) {
    
    // update age_manager if necessary
    age_manager = WeibullAgeSampler( age_shape, age_scale, 100 );
    
    nodes.clear();
    maxID = 0;

    while( !death_events.empty() )
        death_events.pop();
    
    for ( int i = 0; i < N; ++i ) {
        // sample age and lifetime

        pair<int, int> age_data = age_manager.sample_joint_age_residual_eq(*mt);
        int age = age_data.first;

        Individual new_individual = Individual( maxID, i, -age );

        vector<int> t_infect = {};
        t_infect.reserve( 4 );
        
        for ( int nd = 0; nd < 4; ++nd ) {

            int t = -1;
            double s = -log( 1. - uniGen( *mt ) ); // !! or 1 0 u?

            if ( s < sero_lam1 * sero_a0  ) {
                t = int( s / sero_lam1 );
            }
            else {
                t = int( sero_a0 + ( s - sero_lam1 * sero_a0 ) / sero_lam2 );
            }
            
            if ( t < age )
                t_infect.push_back( t - age ); // get timing of infection
            
        }
        
        int n_infect = static_cast<int>( t_infect.size() );
        
        
        if ( n_infect > 0 ) {
            
            // sort infection times
            sort( t_infect.begin(), t_infect.end() );
            
            // get random strains
            vector<int> strains = { 0, 1, 2, 3 };
            std::shuffle( strains.begin(), strains.end(), *mt );
            
            for ( int j = 0; j < n_infect; ++j ) {
                
                int t = t_infect[j];
                int strain = strains[j];
                
                if ( j < n_infect - 1 ) {
                    new_individual.add_immunity( strain, t, true ); // b4 to last infection events
                    new_individual.add_infection_ev( strain, t );
                }
                
                else { // last infection event
                    
                    // check latency
                    int dtEI = compM.draw_latency_time( strain ) ;
                    if ( t + dtEI > 0 ) {
                        
                        new_individual.set_latent_state( true );
                        latency_events.push( InfEvent( new_individual.get_index(), new_individual.get_id(), strain, t + dtEI ) );
                        
                    }
                    else {
                        
                        int dtIR = 1 + floor( log( 1. - uniGen( *mt ) ) / log( 1. - compM.get_rec_prob( strain ) ) ) ;
                        
                        if ( t + dtEI + dtIR > 0 ) {
                            new_individual.apply_infection( t, strain, true );
                            infected_nodes.insert( new_individual.get_index() );
                            
                        }
                        else {
                            int dtRS = compM.draw_waning_time( strain );
                            if ( t + dtEI + dtIR + dtRS > 0 ) {
                                new_individual.add_immunity( strain, t, true );
                                new_individual.add_infection_ev( strain, t );
                                new_individual.update_crossImmD( 1 );
                                
                                waning_events.push( WaneEvent( new_individual.get_index(), new_individual.get_id(), strain, t + dtEI + dtIR + dtRS ) );
                                
                            }
                            else {
                                new_individual.add_immunity( strain, t, true );
                                new_individual.add_infection_ev( strain, t );
                            }
                        
                        }
                        
                    }
                    
                }
                
            }
    
        }
        
        nodes.push_back( new_individual );
        death_events.push( DeathEvent( i, age_data.second ) );
        maxID++;
        
    }
    
}



void Community::do_intro_step(int t, vector<int>& strains_intro, CompartmentalModel &compM) {
    for (int& strain: strains_intro) {
        int n_intros = poi_intros(*mt);
        if (n_intros > 0)
            set_random_infection(strain, n_intros, t, compM);
    }
}

void Community::do_demography_step(int t, CompartmentalModel& compM, bool check_v_tr, bool record_inf) {
    //int count = 0;
    
    //std::geometric_distribution<int> geomtest(1./(20 * 365));
    
    while(death_events.top().t == t) {
        int idx = death_events.top().index;
        if ( nodes[idx].is_infectious() )
            infected_nodes.erase( idx );
        maxID++;
        nodes[idx].death(t, maxID);
        
        death_events.pop();
        death_events.push( DeathEvent( idx, t + age_manager.sample_lifetime(*mt) ) );
        
        // vertical transmission Zika
        if ( check_v_tr ) {
            // draw parent
            // must be sexually active ( age > 15 )
            int idx_parent = 0;
            do
                idx_parent = randint_population(*mt);
            while( nodes[idx_parent].get_age_y(t) >= 15 and nodes[idx_parent].get_age_y(t) < 50 );
            
            auto& target = nodes[idx];
            // check if parent infectious with ZIKV
            if ( nodes[idx_parent].is_infectious_with( 4, compM ) ) {

                if ( record_inf ) // record infection event
                    inf_events_output.emplace_back( target.get_id(), t, target.get_age_y(t), 4, 0, 0, -1, -1, -1 );

                target.apply_infection( t, 4, false );
                infected_nodes.insert( idx );
                
                //std::cout << t << " newborn" << std::endl;
                
            }
            
        }
        //count++;
    }
    //cout << t << " " << count << endl;
}


void Community::check_waning(int t, CompartmentalModel& compM) {
    
    if (waning_events.size() == 0)
        return;
    
    while( waning_events.top().t == t ) {

        int idx = waning_events.top().index;
        int ID = waning_events.top().ID;
        
        if ( nodes[idx].get_id() == ID ) // check ID: individual might have been replaced by newborn
            nodes[idx].update_crossImm( waning_events.top().strain, -1, compM );
        waning_events.pop();
        
        if (waning_events.size() == 0)
            return;
    }
}

void Community::check_latency(int t, CompartmentalModel &compM) {
    if ( latency_events.size() == 0 )
        return;
    
    while ( latency_events.top().t == t ) {
        int idx = latency_events.top().index;
        int ID = latency_events.top().ID;
        int strain = latency_events.top().strain;
        
        if ( nodes[idx].get_id() == ID ) {// check ID: individual might have been replaced by newborn
            if ( nodes[idx].apply_infection( t, strain, compM.is_endemic(strain) ) ) {
                infected_nodes.insert( idx );
            }
            //nodes[idx].update_crossImm( waning_events.top().strain, -1, compM );
        }
        latency_events.pop();
        
        if ( latency_events.size() == 0)
            return;
    }
}

void Community::do_dynamics_step(int t, CompartmentalModel& compM) {
    int deg;
    int idx2;
    double seas_factor;

    //unordered_set<int> to_update (N);
    
    // random sort strains
    compM.shuffle_strain_set();
    
    compute_seasonal_forcing_t(t);
    
    for (const int& idx: infected_nodes) {
        deg = poi_contacts(*mt);
        
        Individual& source = nodes[idx];
        vector<int>& strains_source = compM.get_carried_strains(source);
        
        for (int i = 0; i < deg; ++i) {
            // select individual
            do {
                idx2 = randint_population( *mt );
            }
            while ( idx2 == idx );
            
            Individual& target = nodes[idx2];
            
            for (int& strain: strains_source) {
                
                if ( compM.is_endemic( strain ) )
                    seas_factor = seas_factor_D;
                else
                    seas_factor = seas_factor_Z;
                
                double p_inf = compM.get_immunity( target, strain ) * seas_factor * compM.get_tr(source, strain);

                
                if ( p_inf > 0. ) {
                    if ( target.update_hazard( hazardify_prob( p_inf ) ) ) {
                        
                        compM.update_incidence( strain );

                        latency_events.push( InfEvent( idx2, target.get_id(), strain, t + compM.draw_latency_time( strain ) ) );
                        
                        target.set_latent_state( true );
                        
                    }
                }
            }
        }
        
        // check recovery
        for (int& strain: strains_source) {
            if ( uniGen( *mt ) < compM.get_rec_prob( strain ) ) {
                rec_evs.push( make_pair( idx, strain ) ); // add recovery event
            }
        }
    }
        
    // finally update those who got infected and those who recovered
    // add/remove them to infectious list
    
    /*
    for (auto ev = inf_evs.begin(), ev_end = inf_evs.end(); ev < ev_end; ++ev) {
        if ( nodes[ev->first].apply_infection(t, ev->second, compM.is_endemic( ev->second ) ) )
            infected_nodes.insert(ev->first);
        
        compM.update_incidence(ev->second);
    }
    inf_evs.clear();
     */
    
    // check latency
    check_latency( t, compM );
    
    for (auto ev = rec_evs.begin(), ev_end = rec_evs.end(); ev < ev_end; ++ev) {
        if ( nodes[ev->first].apply_recovery(ev->second) )
            infected_nodes.erase(ev->first);
        nodes[ev->first].update_crossImm(ev->second, 1, compM);
        int delta_t = compM.draw_waning_time(ev->second);
        waning_events.push( WaneEvent(ev->first, nodes[ev->first].get_id(), ev->second, t + delta_t ) );
    }
    rec_evs.clear();
    
    //std::cout << infected_nodes.size() << std::endl;

    
    // process waning of cross-immunity
    check_waning( t, compM );
       
}


void Community::do_dynamics_step_output(int t, CompartmentalModel& compM) {
    int deg;
    int idx2;
    double seas_factor;
    
    //unordered_set<int> to_update (N);
    
    // random sort strains
    compM.shuffle_strain_set();
    
    compute_seasonal_forcing_t(t);

    for (const int& idx: infected_nodes) {
        deg = poi_contacts(*mt);
        
        Individual& source = nodes[idx];
        vector<int>& strains_source = compM.get_carried_strains(source);
        
        for (int i = 0; i < deg; ++i) {
            // select individual
            do {
                idx2 = randint_population(*mt);
            }
            while (idx2 == idx);
            
            Individual& target = nodes[idx2];
            for (int& strain: strains_source) {                
                
                if ( compM.is_endemic( strain ) )
                    seas_factor = seas_factor_D;
                else
                    seas_factor = seas_factor_Z;

                double p_inf = compM.get_immunity(target, strain) * seas_factor * compM.get_tr(source, strain);
                
                if ( p_inf > 0. ) {
                    if ( target.update_hazard( hazardify_prob(p_inf) ) ) {                        
                        //inf_evs.push( make_pair(idx2, strain) );
                        
                        compM.update_incidence( strain );
                        
                        int nExposuresD = target.is_exposed_toD(compM); // count previous exposures to D
                        int nExposuresZ = target.is_exposed_toZ(compM); // count previous exposures to Z
                        
                        // store data
                        inf_events_output.emplace_back(target.get_id(), t, target.get_age_y(t), strain, nExposuresD, nExposuresZ, target.get_last_tD(), target.get_first_tD(), target.get_tZ(), target.get_inf_evs() );
                        
                        latency_events.push( InfEvent( idx2, target.get_id(), strain, t + compM.draw_latency_time( strain ) ) );
                        
                        target.set_latent_state(true);
                    }
                }
                
            }
        }
        
        // check recovery
        for (int& strain: strains_source) {
            if ( uniGen(*mt) < compM.get_rec_prob(strain) ) {
                rec_evs.push( make_pair(idx, strain) ); // add recovery event
            }
        }
    }
        
    // finally update those who got infected and those who recovered
    // add/remove them to infectious list
    
    /*
    for (auto ev = inf_evs.begin(), ev_end = inf_evs.end(); ev < ev_end; ++ev) {
        
        int nExposuresD = nodes[ev->first].is_exposed_toD(compM); // count previous exposures to D
        int nExposuresZ = nodes[ev->first].is_exposed_toZ(compM); // count previous exposures to Z
        
        // store data
        auto& node = nodes[ev->first];
        inf_events_output.emplace_back(node.get_id(), t, node.get_age_y(t), ev->second, nExposuresD, nExposuresZ, node.get_last_tD(), node.get_first_tD(), node.get_tZ() );
        
        
        if ( nodes[ev->first].apply_infection(t, ev->second,  compM.is_endemic( ev->second ) ) )
            infected_nodes.insert(ev->first);
        
        compM.update_incidence(ev->second);
        
    }
    inf_evs.clear();
    */
    
    // check_latency
    check_latency( t, compM );
    
    for (auto ev = rec_evs.begin(), ev_end = rec_evs.end(); ev < ev_end; ++ev) {
        if ( nodes[ev->first].apply_recovery(ev->second) )
            infected_nodes.erase(ev->first);
        
        nodes[ev->first].update_crossImm( ev->second, 1, compM );
        int delta_t = compM.draw_waning_time( ev->second );
        waning_events.push( WaneEvent(ev->first, nodes[ev->first].get_id(), ev->second, t + delta_t ) );
    }
    rec_evs.clear();
    
    // process waning of cross-immunity
    check_waning(t, compM);
}


void Community::initialize_immunity_strain(int strain, CompartmentalModel &compM) {
    double R0 = compM.get_tr(strain) * mean_degree / compM.get_rec_prob(strain);
    
    if (R0 <= 1.)
        return;
    
    double p = 1 - 1./R0; // probability of being exposed
    
    geometric_distribution<int> geom( p );
    
    int n = 1 + geom(*mt);
    
    while (n <= N) {
        nodes[n - 1].add_immunity(strain);
        n += 1 + geom(*mt);
    }
}

void Community::set_random_infection(int strain, int n, int t, CompartmentalModel& compM) {
    int idx = 0;
    int n_attempts = 0;
    bool is_endemic = compM.is_endemic( strain );
    for (int i = 0; i < n; ++i) {
        idx = randint_population(*mt);
        if ( (!nodes[idx].is_immune_to(strain, compM) ) and ( !nodes[idx].is_latent() ) ) {
            n_attempts = 0;
            if ( nodes[idx].apply_infection(t, strain, is_endemic ) )
                infected_nodes.insert(idx);
        }
        else {
            n_attempts++;
            i--;
        }
        
        if (n_attempts == MAX_INTRO_ATTEMPTS)
            break;
    }
}

void Community::compute_seasonal_forcing_t(const int& t) {
    
    seas_factor_D = 1 + seasonal_force_D * cos( ( t - seasonal_max_t_D ) * 2 * PI/365. );
    seas_factor_Z = 1 + seasonal_force_Z * cos( ( t - seasonal_max_t_Z ) * 2 * PI/365. );

}

vector<int> Community::get_prevalence_strains(CompartmentalModel &compM) {
    vector<int> histo_state = vector<int>( (int)pow(2, compM.get_nstrains() ), 0);
    vector<int> histo_strain = vector<int>( compM.get_nstrains() , 0);
    
    for (int idx: infected_nodes) {
        Individual& indiv = nodes[idx];
        histo_state[ indiv.get_inf_state() ]++;
    }
    
    int state = 0;
    for (int count: histo_state) {
        if ( count > 0 ) {
            for (int strain: compM.get_carried_strains(state) ) {
                histo_strain[strain] += count;
            }
        }
        state++;
    }
    
    return histo_strain;
}

// compute age profile histogram
vector<int> Community::get_age_structure(int t) {
    vector<int> histo = vector<int>(101, 0);
    int age;
    for (Individual& ind: nodes) {
        age = (int) ( ( t - ind.get_birth() ) / 365 );
        if (age >= 100)
            histo[101]++;
        else
            histo[age]++;
    }    
    return histo;
}

void Community::flush_events(ostream& outfile) {
    for (auto& ev: inf_events_output) {
        outfile << ev.id << "," << ev.t << "," << ev.age  << "," << ev.strain << endl;
    }
    inf_events_output.clear();
}

void Community::print_AR_age_endemic(int t, CompartmentalModel &compM, ostream &outfile) {
    int nmax = compM.get_nstrains() - 1;
    vector<int> histo = vector<int>(101, 0);
    for (auto& ev: inf_events_output) {
        if (ev.strain < nmax) // select endemic strain infections only
            histo[ev.age]++;
    }
    outfile << t << "," << stringify_vector(histo) << endl;
}

// prints 3 histograms: one for primary exposures, one for secondary exposures and one for other types
void Community::print_AR_age_endemic_exp(int t, CompartmentalModel &compM, ostream &outfile) {
    int nmax = compM.get_nstrains() - 1;
    
    // select endemic strain infections only
    
    // select primary exposures only
    {
        vector<int> histo = vector<int>(101, 0);
        for (auto& ev: inf_events_output) {
            if ( (ev.strain < nmax) and (ev.previousExpD == 0) )
                histo[ev.age]++;
        }
        outfile << t << ",0," << stringify_vector(histo) << endl;
    }
    
    // select secondary exposures only
    {
        vector<int> histo = vector<int>(101, 0);
        for (auto& ev: inf_events_output) {
            if ( (ev.strain < nmax) and (ev.previousExpD == 1) )
                histo[ev.age]++;
        }
        outfile << t << ",1," << stringify_vector(histo) << endl;
    }
    
    // select post-secondary exposures only
    {
        vector<int> histo = vector<int>(101, 0);
        for (auto& ev: inf_events_output) {
            if ( (ev.strain < nmax) and (ev.previousExpD > 1) )
                histo[ev.age]++;
        }
        outfile << t << ",2," << stringify_vector(histo) << endl;
    }
}

void Community::print_AR_age_invasive_exp(int t, CompartmentalModel &compM, ostream &outfile) {
    int nmax = compM.get_nstrains() - 1;
    
    // select primary exposures only
    {
        vector<int> histo = vector<int>(101, 0);
        for (auto& ev: inf_events_output) {
            if ( (ev.strain == nmax) and (ev.previousExpD == 0) )
                histo[ev.age]++;
        }
        outfile << t << ",0," << stringify_vector(histo) << endl;
    }
    
    // select secondary exposures only
    {
        vector<int> histo = vector<int>(101, 0);
        for (auto& ev: inf_events_output) {
            if ( (ev.strain == nmax) and (ev.previousExpD == 1) )
                histo[ev.age]++;
        }
        outfile << t << ",1," << stringify_vector(histo) << endl;
    }
    
    // select post-secondary exposures only
    {
        vector<int> histo = vector<int>(101, 0);
        for (auto& ev: inf_events_output) {
            if ( (ev.strain == nmax) and (ev.previousExpD > 1) )
                histo[ev.age]++;
        }
        outfile << t << ",2," << stringify_vector(histo) << endl;
    }
    
    
}

void Community::print_AR_age_invasive(int t, CompartmentalModel &compM, ostream &outfile) {
    int nmax = compM.get_nstrains() - 1;
    vector<int> histo = vector<int>(101, 0);
    for (auto& ev: inf_events_output) {
        if (ev.strain == nmax) // select invasive strain infections only
            histo[ev.age]++;
    }
    outfile << t << "," << stringify_vector(histo) << endl;
}

// prints 6 histograms, splitting cases according to 0,1,>1 past exposures to Dengue and exposure to Z
void Community::print_AR_age_endemic_expZ(int t, CompartmentalModel &compM, ostream &outfile) {
    int nmax = compM.get_nstrains() - 1;
    
    // help
    int ct = 0;
    for (auto& ev: inf_events_output) {
        if ( ev.previousExpZ == 0 )
            ct++;
    }
    
    // n.b. select endemic strain infections only
    
    // A) No exposure to Zika histograms
    
    // select primary exposures only
    {
        vector<int> histo = vector<int>(101, 0);
        for (auto& ev: inf_events_output) {
            if ( (ev.strain < nmax) and (ev.previousExpD == 0) and (ev.previousExpZ == 0) )
                histo[ev.age]++;
        }
        outfile << t << ",0,0," << stringify_vector(histo) << endl;
    }
    
    // select secondary exposures only
    {
        vector<int> histo = vector<int>(101, 0);
        for (auto& ev: inf_events_output) {
            if ( (ev.strain < nmax) and (ev.previousExpD == 1) and (ev.previousExpZ == 0) )
                histo[ev.age]++;
        }
        outfile << t << ",0,1," << stringify_vector(histo) << endl;
    }
    
    // select post-secondary exposures only
    {
        vector<int> histo = vector<int>(101, 0);
        for (auto& ev: inf_events_output) {
            if ( (ev.strain < nmax) and (ev.previousExpD > 1) and (ev.previousExpZ == 0) )
                histo[ev.age]++;
        }
        outfile << t << ",0,2," << stringify_vector(histo) << endl;
    }
    
    // B) With exposure to Zika histograms

    // select primary exposures only
    {
        vector<int> histo = vector<int>(101, 0);
        for (auto& ev: inf_events_output) {
            if ( (ev.strain < nmax) and (ev.previousExpD == 0) and (ev.previousExpZ > 0) )
                histo[ev.age]++;
        }
        outfile << t << ",1,0," << stringify_vector(histo) << endl;
    }
    
    // select secondary exposures only
    {
        vector<int> histo = vector<int>(101, 0);
        for (auto& ev: inf_events_output) {
            if ( (ev.strain < nmax) and (ev.previousExpD == 1) and (ev.previousExpZ > 0) )
                histo[ev.age]++;
        }
        outfile << t << ",1,1," << stringify_vector(histo) << endl;
    }
    
    // select post-secondary exposures only
    {
        vector<int> histo = vector<int>(101, 0);
        for (auto& ev: inf_events_output) {
            if ( (ev.strain < nmax) and (ev.previousExpD > 1) and (ev.previousExpZ > 0) )
                histo[ev.age]++;
        }
        outfile << t << ",1,2," << stringify_vector(histo) << endl;
    }
    
}


void Community::print_age_distribution(int t, string outfile_path) {
    vector<int> histo = vector<int>(101, 0);
    int a;
    for (Individual& ind: nodes) {
        a = min( static_cast<int>( ( t - ind.get_birth() ) / 365. ), 100 );
        histo[a]++;
    }
    ofstream outfile;
    outfile.open(outfile_path);
    outfile << stringify_vector(histo);
    outfile.close();
}

void Community::print_seroprevalence(int t, int n, bool checkDengue, bool checkZika, CompartmentalModel& compM, ostream& outfile) {
    
    if ( !(checkDengue or checkZika) )
        return; // no check to be made
    
    string SeroCheck = "";
    if (checkDengue)
        SeroCheck += "D";
    if (checkZika)
        SeroCheck += "Z";
    
    vector<int> histo = vector<int>(101, 0);
    
    for (int i = 0; i < n; ++i) {
        int idx = randint_population(*mt);
        bool hasD = false;
        bool hasZ = false;
        for ( int strain: compM.get_imm_strains(nodes[idx]) ) {
            if ( compM.is_invasive(strain) )
                hasZ = true;
            else
                hasD = true;
        }
        
        if (checkDengue) {
            if (!hasD)
                continue;
        }
        
        if (checkZika) {
            if (!hasZ)
                continue;
        }
        
        int a = min( static_cast<int>( ( t - nodes[idx].get_birth() ) / 365. ), 100 );
        histo[a]++;
        
    }
    
}

// get seroprevalence for Dengue
// create separate histograms according to the number of strains seen
void Community::print_seroprevalenceD(int t, CompartmentalModel& compM, ostream& outfile) {
    int nstrainsD = compM.get_nstrains() - 1;
    
    std::vector<std::vector<int>> histo_vec;
    for (int i = 0; i < nstrainsD; i++) {
        histo_vec.push_back( std::vector<int>(101, 0) );
    }
    
    for (Individual& ind: nodes) {
        int count = -1;
        for ( int& strain: compM.get_imm_strains(ind) ) {
            if ( compM.is_endemic(strain) )
                count++;
        }
        
        if ( count >= 0 ) {
            int a = min( static_cast<int>( ( t - ind.get_birth() ) / 365. ), 100 );
            histo_vec[count][a]++;
        }
    }
    
    for (int i = 1; i <= nstrainsD; i++) {
        outfile << t << "," << i << "," << stringify_vector(histo_vec[i - 1]) << endl;
    }
}

void Community::check_prevalence() {
    int count = 0;
    for (Individual& ind: nodes) {
        if ( ind.is_infectious() )
            count++;
    }
    assert( count == infected_nodes.size() );
}


// loop over individuals and print time since last Dengue infection and age (in days)
void Community::print_TlastDengueInfection(int t, CompartmentalModel &compM, string outfile_path) {
    
    ofstream outfile;
    outfile.open(outfile_path);
    
    bool hadD;
    int a;
    int t_inf;
    int strain;
    
    for (Individual& ind: nodes) {
        hadD = false;
        auto& evs = ind.get_inf_evs();
        
        if ( evs.size() == 0 ) {
            a = t - ind.get_birth();
            outfile << a << "," << -1 << std::endl;
            continue;
        }

        for (auto it = evs.rbegin(); it != evs.rend(); ++it) {
            t_inf = it->t;
            strain = it->strain;
            if ( compM.is_endemic( strain ) ) {
                a = t - ind.get_birth();
                outfile << a << "," << t - t_inf << std::endl;
                hadD = true;
                break;
            }
        }
        
        if (!hadD) {
            a = t - ind.get_birth();
            outfile << a << "," << -1 << std::endl;
        }
    }
    
    outfile.close();
}

void Community::print_infection_events(int t, CompartmentalModel &compM, std::string outfile_path) {

    ofstream outfile;
    //outfile.open(outfile_path, std::ofstream::out | std::ofstream::app);
    outfile.open(outfile_path);

    outfile << "t,strain,age,nD,nZ,tD,t1D,tZ" << std::endl;
    for (auto& ev: inf_events_output) {
        outfile << ev.t << "," << ev.strain << "," << ev.age << "," << ev.previousExpD << "," << ev.previousExpZ << "," << ev.tLastD << "," << ev.tFirstD << "," << ev.tZ << std::endl;
    }

    outfile.close();
}

void Community::print_infection_events_flush(int t, CompartmentalModel &compM, ostream& outfile) {
    for (auto& ev: inf_events_output) {
        outfile << ev.t << "," << ev.strain << "," << ev.age << "," << ev.previousExpD << "," << ev.previousExpZ << "," << ev.tLastD << "," << ev.tFirstD << "," << ev.tZ << std::endl;
    }
    inf_events_output.clear();
}


void Community::check_naive_Z(CompartmentalModel& compM) {
    uint ct = 0;
    for (auto& indiv: nodes) {
        if ( indiv.is_exposed_toZ( compM ) == 0 ) {
            ct++;
        }
    }
    cout << "# non-exposed to Z: " << ct << endl;
}

void Community::fill_agehisto_secondaryDengue( vector<int>& histo, CompartmentalModel &compM ) {
    int amax = static_cast<int>( histo.size() ) - 1; // censor age to amax
    int age = 0;
    for ( auto& ev: inf_events_output ) {
        if ( compM.is_endemic( ev.strain ) ) {
            age = ( ev.age <= amax ) ? ev.age : amax;
            if ( ev.previousExpD == 1 ) {
                ++histo[ age ];
            }
        }
    }
}

void Community::fill_agehisto_parity( vector<vector<int>>& histo, CompartmentalModel &compM ) {
    int amax = static_cast<int>( histo[0].size() ) - 1; // censor age to amax
    int age = 0;
    int parity = 0;
    for ( auto& ev: inf_events_output ) {
        if ( compM.is_endemic( ev.strain ) ) {
            parity = ev.previousExpD;
            age = ( ev.age <= amax ) ? ev.age : amax;
            ++histo[ parity ][ age ];
        }
    }
}

void Community::fill_agehisto_parity_singleDvsAll( vector<vector<vector<int>>>& histo, CompartmentalModel& compM )  {
    int amax = static_cast<int>( histo[0].size() ) - 1; // censor age to amax
    int age = 0;
    int parity = 0;
    for ( auto& ev: inf_events_output ) {
        if ( compM.is_endemic( ev.strain ) ) {
            parity = ev.previousExpD;
            age = ( ev.age <= amax ) ? ev.age : amax;
            if ( ev.strain == 0 )
                ++histo[ parity ][ age ][ 0 ]; // single strain only
            ++histo[ parity ][ age ][ 1 ]; // all endemic strains together
        }
    }
}

void Community::fill_agehisto_parity_Z( vector<vector<int>>& histo, CompartmentalModel &compM ) {
    int amax = static_cast<int>( histo[0].size() ) - 1; // censor age to amax
    int age = 0;
    int parity = 0;
    
    int idx = 0;
    
    for ( auto& ev: inf_events_output ) {
        if ( compM.is_invasive( ev.strain ) ) {
            parity = ev.previousExpD;
            
            if ( parity == 0 )
                idx = 0;
            else if ( parity == 1 )
                idx = 1;
            else
                idx = 2;
            
            age = ( ev.age <= amax ) ? ev.age : amax;
            ++histo[ idx ][ age ];
        }
    }
}


void Community::fill_agehisto_parity_ZD( vector<vector<int>>& histo, CompartmentalModel &compM ) {
    int amax = static_cast<int>( histo[0].size() ) - 1; // censor age to amax
    int age = 0;
    int parityD = 0;
    int parityZ = 0;

    
    int idx = 0;
    
    for ( auto& ev: inf_events_output ) {
        // Z infection
        if ( compM.is_invasive( ev.strain ) ) {
            parityD = ev.previousExpD;
            
            if ( parityD == 0 )
                idx = 0;
            else if ( parityD == 1 )
                idx = 1;
            else if ( parityD == 2 )
                idx = 2;
            else
                idx = 3;
        }
        else { // D infection
            parityD = ev.previousExpD;
            parityZ = ev.previousExpZ;
            
            if ( parityD == 0 )
                idx = 4;
            else if ( parityD == 1 )
                idx = 5;
            else if ( parityD == 2 )
                idx = 6;
            else
                idx = 7;
            
            if ( parityZ == 1 )
                idx += 4;
        }
        
        age = ( ev.age <= amax ) ? ev.age : amax;
        ++histo[ idx ][ age ];
    }
}


// counts the number of seropositive individuals by age and parity of Dengue exposure status
void Community::fill_seroprevalenceD( int t, vector<vector<int>>& histo, CompartmentalModel& compM ) {
    
    for ( Individual& ind: nodes ) {
    
        int nD = ind.get_nD();
        int a  = min( static_cast<int>( ( t - ind.get_birth() ) / 365. ), 100 );
        histo[nD][a]++;
    
    }
    
}

