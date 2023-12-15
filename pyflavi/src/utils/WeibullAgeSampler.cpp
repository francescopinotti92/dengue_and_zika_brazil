//
//  WeibullAgeSampler.cpp
//  CrossFlavivirus
//
//  Created by MacBook Pro on 07/01/2021.
//

#include "WeibullAgeSampler.hpp"


WeibullAgeSampler::WeibullAgeSampler(double shape_, double scale_, int max_age_): shape(shape_), scale(scale_ * 365.), max_age(max_age_ * 365.) {
    
    _uniGen = std::uniform_real_distribution<double>(0., 1.);
    
}


// sample lifetime for a newborn (in days)
int WeibullAgeSampler::sample_lifetime(std::mt19937_64 &mt_) {
    double u = 0.;
    int l = 0;
    do {
        u = _uniGen(mt_);
        l = static_cast<int>( scale * pow( -log( 1 - u ), 1./shape ) );
    }
    while ( l >= max_age );
   
    return l;
}

int WeibullAgeSampler::sample_residual_eq(std::mt19937_64 &mt_) {
    
    // use acceptance-rejection
    double u = 0.;
    double l_proposal = 0;
    double p_accept = 0.;
    
    // compute acceptance-rejection constant
    double M = 0.;
    if ( shape == 1 )
        M = 1.;
    else {
        M = exp( ( 1 - 1./shape ) * pow( 1/shape, 1 / ( shape - 1 ) ) );
    }
    
    while ( true ) {
        // sample from proposal distribution
        do {
            u = _uniGen(mt_);
            l_proposal = -log( 1 - u ) * scale;
        }
        while ( l_proposal >= max_age );
        
        // accept proposal?
        p_accept = exp( l_proposal / scale - pow( l_proposal / scale, shape ) ) / M;
        if ( _uniGen(mt_) < p_accept )
            return static_cast<int>( l_proposal );
    }
}

// sample age and residual lifetime distribution to initialize the system (both in days)
std::pair<int, int> WeibullAgeSampler::sample_joint_age_residual_eq(std::mt19937_64 &mt_) {
    int age = 0;
    int residual = 0;
    int lifetime = 0;
    
    // sample age from marginal distribution
    age = sample_residual_eq( mt_ ); // eq age and residual have same marginal distributions
    
    // sample lifetime, conditioned on being >= than age
    if ( age == max_age ) {
        lifetime = age;
    }
    else {
        do {
            lifetime = sample_lifetime( mt_ );
        }
        while ( lifetime < age );
    }
   
    residual = lifetime - age;

    
    return std::make_pair( age, residual );
}
