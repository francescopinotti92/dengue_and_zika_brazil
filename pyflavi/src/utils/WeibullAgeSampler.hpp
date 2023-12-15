//
//  WeibullAgeSampler.hpp
//  CrossFlavivirus
//
//  Created by MacBook Pro on 07/01/2021.
//

#ifndef WeibullAgeSampler_hpp
#define WeibullAgeSampler_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <random>

class WeibullAgeSampler {
public:
    WeibullAgeSampler(double shape, double scale, int max_age); // scale and max_age should be given in years
    int sample_lifetime(std::mt19937_64& mt_);
    int sample_residual_eq(std::mt19937_64& mt_);
    std::pair<int, int> sample_joint_age_residual_eq(std::mt19937_64& mt_);
private:
    double shape;
    double scale;
    int max_age;
    std::uniform_real_distribution<double> _uniGen;
};



struct DeathEvent {
    DeathEvent(int index_, int t_): index(index_), t(t_) {};
    int index;
    int t;
};

struct CompareDeathEvents {
    CompareDeathEvents() {};
    inline bool operator() (const DeathEvent& left, const DeathEvent& right) {
        return left.t > right.t;
    }
};

#endif /* WeibullAgeSampler_hpp */
