//
//  individual.hpp
//  CrossFlavivirus
//
//  Created by MacBook Pro on 05/01/2021.
//

#ifndef individual_hpp
#define individual_hpp

#define MAXAGE_DEATHPROB 90.0
#define SHAPE_DEATHPROB 12.0
#define PI 3.14159265
#define MAX_INTRO_ATTEMPTS 5000

#include "event_vec.hpp"
#include "WeibullAgeSampler.hpp"
#include <algorithm>
#include <limits>
#include <cassert>
#include <cmath>
#include <stdio.h>
#include <fstream>
#include <bitset>
#include <random>
#include <string>
#include <queue>
#include <vector>
#include <unordered_set>

using namespace std;

extern mt19937_64* mt;
extern uniform_int_distribution<int> randint_dayofyear;
extern uniform_int_distribution<int> randint_population;
extern uniform_real_distribution<double> uniGen;
extern exponential_distribution<double> expoGen;
//geometric_distribution<int> geoGen;

//extern uniform_int_distribution<int> randint_dayofyear(0, 364);
//extern uniform_int_distribution<int> randint_population(0, 1);

void assign_mt(mt19937_64* mt_);

template <typename T>
string stringify_vector(const vector<T> &v, string last_delim = "") {
    string res;
    string delim;
    if (v.size() > 0) {
        for (unsigned int ii = 0; ii < v.size(); ++ii) {
            delim = (ii < v.size() - 1) ? "," : last_delim;
            res += to_string(v[ii]) + delim;
        }
    }
    else
        res = last_delim;
    return res;
}

inline double hazardify_prob(const double & p) {
    return -std::log(1. - p);
}

struct AgeSamplerClass {
    AgeSamplerClass(int age_max_ = 100);
    vector<double> psurv;
    vector<double> psurv_cum;
    int age_max;
};

double ageDependentMortality( const double& age);
int sampleAge();

struct IndividualInfectionEvent {
    IndividualInfectionEvent(): strain( -1 ), t( std::numeric_limits<int>::min() ) {};
    IndividualInfectionEvent( const int& strain, const int& t ): strain( strain ), t( t ) {};
    int strain;
    int t;
};

class Individual;

class CompartmentalModel {
public:
    CompartmentalModel(int nstrains_, vector<double> tr, vector<double> sig, vector<double> ltcy );
        
    void set_immunity_params( double lD, double lZ, double rhoDZ, double rhoZD );
    void set_immunity_params( double lD, double lZ, double rhoDZ, double rhoZD, double rhoZD0, double csiDZ, double csiZD, double gamma1D, double gamma1Z, double gamma2D, double gamma2Z );

    void set_immunity_params_v2( double lD, double lZ, double rhoZD0, double rhoZD1, double rhoZD2, double chiZD0, double chiZD1, double chiZD2 );
    
    inline int    draw_latency_time( const int& strain ) { return 1 + latencyGens[strain](*mt); }
    inline int    draw_waning_time( const int& strain ) { return 1 + waningGens[strain](*mt); }
    inline int    get_nstrains() { return nstrains; }
    inline int get_Z_id() { return nstrains - 1; }
    inline double get_tr(const int& strain) { return tr_i[strain]; }
    inline void set_tr( const int& strain, const double& value ) { tr_i[strain] = value; }
    double get_tr( Individual& source, const int& strain ) ;
    inline double get_sigma(const int& strain) { return sig_i[strain]; }
    inline double get_rec_prob(const int& strain) { return rec_prob_i[strain]; }
    inline bool is_invasive(const int& strain) { return ( (strain == nstrains - 1) ? true : false); }
    inline bool is_endemic(const int& strain) { return !is_invasive(strain); }
    inline bool is_immune( const int& strain, const Individual& indiv );
    //double get_hazard(Individual& target, int strain);
    double get_immunity(Individual &target, const int& strain);
    inline vector<int>& get_carried_strains(Individual& indiv);
    inline vector<int>& get_carried_strains(const int& inf_state);

    inline vector<int>& get_imm_strains(Individual& indiv);
    
    void shuffle_strain_set();
    
    inline void update_incidence(const int& strain) { cum_incidence[strain]++; }
    void clear_incidence() {
        for (int i = 0; i < nstrains; i++)
            cum_incidence[i] = 0;
    }
    void print_cumulative_incidence(int t, ostream& outfile);
    
    
private:
    int nstrains;
    
    vector<double> tr_i;
    vector<double> sig_i;
    vector<double> rec_prob_i;
    vector<double> ltc_prob_i;
    
    /*
    double rhoDZ;       // short-term D-induced mod. susc. against Z
    double rhoZD;       // short-term Z-induced mod. susc. against post-primary D
    double rhoZD0;      // short-term Z-induced mod. susc. against primary D
    double gamma1D;     // long-term D1-induced mod. susc. against D
    double gamma1Z;     // long-term D1-induced mod. susc. against Z
    double gamma2D;     // long-term D2+-induced mod. susc. against D
    double gamma2Z;     // long-term D2+-induced mod. susc. against Z
    */
    
    double csiDZ;       // short-term D-induced mod. tr. against Z
    double csiZD;       // short-term Z-induced mod. tr. against Z
    
    vector<double> gamma_noZ_D; // relative susc. to D w/out Z
    vector<double> gamma_ltZ_D; // relative susc. to D w/ Z (transient imm. waned)
    vector<double> gamma_stZ_D; // relative susc. to D w/ Z (transient imm. active)
    vector<double> gamma_ltD_Z; // relative susc. to Z w/ Z (transient imm. waned)
    vector<double> gamma_stD_Z; // relative susc. to Z w/ Z (transient imm. active)
    
    vector<vector<double>> cross_imm_susc;
    vector<vector<int>> strains_carried;
    vector<vector<bool>> strains_carried_bool;
    vector<int> cum_incidence;
    
    vector<geometric_distribution<int>> waningGens;
    vector<geometric_distribution<int>> latencyGens;


};

class Individual {
public:
    Individual( int id = 0, int index = 0, int dayborn = 0 );
    
    
    bool apply_infection(const int& t, const int& strain, const bool& isD);
    bool apply_recovery(const int& strain);
    void add_immunity(const int& strain);
    void add_immunity(const int& strain, const int& t, const bool& isD);
    
    void update_crossImm(const int& strain, const int& amount, CompartmentalModel& compM);

    inline void update_crossImmD(const int& amount) { imm_ctrD += amount; } // should be -1 or +1
    inline void update_crossImmZ(const int& amount) { imm_ctrZ += amount; } // should be -1 or +1
    
    void death(int t, int new_id);
    bool update_hazard(double hz);
    inline int get_id() { return id; }
    inline int get_index() { return index; }
    inline int get_birth() { return dayborn; }
    inline int get_age_y(int t) { return (int)( (t - dayborn)/365 ); }
    inline int get_inf_state() const { return inf_state; }
    inline int get_imm_state() const { return imm_state; }
    inline bool is_infectious() { return (inf_state == 0) ? false : true; }
    inline int get_ctrD() { return imm_ctrD; }
    inline int get_ctrZ() { return imm_ctrZ; }
    inline int get_last_tD() { return tLastD; }
    inline int get_first_tD() { return tFirstD; }
    inline int get_tZ() { return tZ; }
    inline int get_nD() { return nD; }
    inline int get_nZ() { return nZ; }
    inline bool is_latent() { return isL; }
    inline void set_latent_state( bool state ) { isL = state; }

    inline void add_infection_ev( const int& strain, const int& t ) { infection_evs.emplace_back( strain, t ); }
    inline vector<IndividualInfectionEvent>& get_inf_evs() { return infection_evs; }

    bool is_infectious_with(const int& strain, CompartmentalModel& compM);
    bool is_immune_to(const int& strain, CompartmentalModel& compM);
    int is_exposed_toD(CompartmentalModel& compM);
    int is_exposed_toZ(CompartmentalModel& compM);
    

private:
    int id;
    int index;
    int dayborn;
    int inf_state;
    int imm_state;
    int nD;
    int nZ;
    int imm_ctrD;
    int imm_ctrZ;
    int tLastD;
    int tFirstD;
    int tZ;
    double hazard;
    bool was_updated;
    bool isL;
    vector<IndividualInfectionEvent> infection_evs;
};


struct EventInfectionO {
    EventInfectionO( const int& id, const int& t, const int& age, const int& strain, const int& pastExpD, const int& pastExpZ, const int& tLastD, const int& tFirstD, const int& tZ, const std::vector<IndividualInfectionEvent>& past_evs = {} ): id( id ), t( t ), age( age ), strain( strain ), previousExpD( pastExpD ), previousExpZ( pastExpZ ), tLastD( tLastD ), tFirstD( tFirstD ), tZ( tZ ), past_evs( past_evs ) {};    
    int id;
    int t;
    int age;
    int strain;
    int previousExpD;
    int previousExpZ;
    int tLastD;
    int tFirstD;
    int tZ;
    std::vector<IndividualInfectionEvent> past_evs;
};


struct WaneEvent {
    WaneEvent(int index_, int ID_, int strain_, int t_): index(index_), ID(ID_), strain(strain_), t(t_) {};
    int index;
    int ID;
    int strain;
    int t;
};

struct CompareWaneEvents {
    CompareWaneEvents() {};
    inline bool operator() (const WaneEvent& left, const WaneEvent& right) {
        return left.t > right.t;
    }
};

struct InfEvent {
    InfEvent(int index_, int ID_, int strain_, int t_): index(index_), ID(ID_), strain(strain_), t(t_) {};
    int index;
    int ID;
    int strain;
    int t;
};

struct CompareInfEvents {
    CompareInfEvents() {};
    inline bool operator() (const InfEvent& left, const InfEvent& right) {
        return left.t > right.t;
    }
};

class Community {
public:
    Community( int N_, double intro_rate_, double mean_degree_, double seasonal_force_D, double seasonal_force_Z, double seasonal_max_t_D, double seasonal_max_t_Z );
    void initialize_population_ageonly( double age_scale, double age_shape );
    void initialize_population_serology( double age_scale, double age_shape, double sero_lam1, double sero_lam2, double sero_a0, CompartmentalModel& compM );
    
    void do_intro_step(int t, vector<int>& strains_intro, CompartmentalModel& compM);
    void do_demography_step(int t, CompartmentalModel& compM, bool check_v_tr = false, bool record_inf = false);
    void do_dynamics_step(int t, CompartmentalModel& compM);
    void do_dynamics_step_output(int t, CompartmentalModel& compM); // stores infection list
    void check_latency(int t, CompartmentalModel& compM);
    void check_waning(int t, CompartmentalModel& compM);
    void initialize_immunity_strain(int strain, CompartmentalModel& compM);
    void set_random_infection(int strain, int n, int t, CompartmentalModel& compM);

    void compute_seasonal_forcing_t(const int& t);
    
    // output functions
    int get_prevalence() { return static_cast<int>( infected_nodes.size() ); }
    void swap_infection_events( vector<EventInfectionO>* v ) {
        v->swap( this->inf_events_output );
    }
    vector<int> get_prevalence_strains(CompartmentalModel& compM);
    vector<int> get_age_structure(int t);
    void flush_events(ostream& outfile);
    void flush_events() { inf_events_output.clear(); }
    void print_AR_age_endemic(int t, CompartmentalModel& compM, ostream& outfile);
    void print_AR_age_endemic_exp(int t, CompartmentalModel& compM, ostream& outfile);
    void print_AR_age_invasive(int t, CompartmentalModel& compM, ostream& outfile);
    void print_AR_age_invasive_exp(int t, CompartmentalModel& compM, ostream& outfile);
    void print_AR_age_endemic_expZ(int t, CompartmentalModel& compM, ostream& outfile);

    void print_age_distribution(int t, string outfile_path);
    void print_TlastDengueInfection(int t, CompartmentalModel& compM, string outfile_path);
    void print_seroprevalence(int t, int n, bool checkDengue, bool checkZika, CompartmentalModel& compM, ostream& outfile);
    void print_seroprevalenceD(int t, CompartmentalModel& compM, ostream& outfile);
    void print_infection_events(int t, CompartmentalModel &compM, std::string outfile_path);
    void print_infection_events_flush(int t, CompartmentalModel &compM, ostream& outfile);
    
    // output functions used with pyflavi
    void fill_agehisto_secondaryDengue( vector<int>& histo, CompartmentalModel& compM );
    void fill_agehisto_parity( vector<vector<int>>& histo, CompartmentalModel& compM );
    void fill_agehisto_parity_singleDvsAll( vector<vector<vector<int>>>& histo, CompartmentalModel& compM );

    
    void fill_agehisto_parity_Z( vector<vector<int>>& histo, CompartmentalModel& compM );
    void fill_agehisto_parity_ZD( vector<vector<int>>& histo, CompartmentalModel& compM );
    
    void fill_seroprevalenceD( int t, vector<vector<int>>& histo, CompartmentalModel& compM );

    int get_n_events() { return static_cast<int>( inf_events_output.size() ); }
    
    // testing functions
    
    void check_prevalence();
    void check_naive_Z(CompartmentalModel& compM);
    
    
private:
    int N;
    int maxID;
    vector<Individual> nodes;
    poisson_distribution<int> poi_contacts;
    poisson_distribution<int> poi_intros;

    unordered_set<int> infected_nodes;
    event_vec<pair<int, int>> inf_evs;
    event_vec<pair<int, int>> rec_evs;
    WeibullAgeSampler age_manager;
    priority_queue<DeathEvent, vector<DeathEvent>, CompareDeathEvents> death_events;
    priority_queue<WaneEvent, vector<WaneEvent>, CompareWaneEvents> waning_events;
    priority_queue<InfEvent, vector<InfEvent>, CompareInfEvents> latency_events;

    vector<EventInfectionO> inf_events_output;

    double intro_rate;
    double mean_degree;
    double seasonal_force_D;
    double seasonal_force_Z;
    double seasonal_max_t_D;
    double seasonal_max_t_Z;
    double seas_factor_D;
    double seas_factor_Z;
    
    
};





#endif /* individual_hpp */
