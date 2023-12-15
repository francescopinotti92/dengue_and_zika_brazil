//
//  main.cpp
//  CrossFlavivirus
//
//  Created by user on 24/08/2021.
//

#include <iostream>
#include "pyflavi.hpp"

// this is just to test code

int main() {
    /*
    simulate_DENV(5, 0, 365 * 101, 365 * 100, 300000, 5, 1.5 / 4.5, 1. / 15.9, 1. / 4.5, 1., 2., 1./30, 0.3, 90., 7., 66., {0,1,2,3}, {365 * 100 + 364});
    std::cout << "done" << std::endl;
     */
    
    int Nsim = 20;
    int Tremoval = 365 * 100;
    int N = 300000;
    int nstrains = 5;
    double epsi = 1./15.9;
    double sigma = 1./4.5;
    double mean_degree = 1.;
    double beta = 1.5 * sigma / mean_degree;
    double lD = 2.;
    double seasonal_force = 0.3;
    double seasonal_peak = 45.;
    double age_shape = 7.;
    double age_max_prob = 66.;
    double intro_rate = 0.1/30;
    
    double lam1 = 0.0001;
    double lam2 = 0.00000392;
    double sero_a0 = 16 * 365;
    
    
    int seed = 0;
    
    simulate_DENV_sero( Nsim, seed, Tmax, Tsample, N, nstrains, beta, epsi, sigma, mean_degree, lD, intro_rate, seasonal_force, seasonal_peak, age_shape, age_max_prob, sero_lam1, sero_lam2, sero_a0, {365 * 100 + 364});
    
    std::cout << "done" << std::endl;
}
