//
//  pyFlavi.hpp
//  CrossFlavivirus
//
//  Created by user on 23/08/2021.
//

#ifndef pyFlavi_h
#define pyFlavi_h

#include <random>
#include <vector>
#include "../../src/model/individual.hpp"


/*
 Simulate the spread of DENV (no ZIKV). Yields the age distribution of Dengue between
 a set of user-defined checkpoints, distinguishing by infection parity
 */
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
                   const std::vector<int>& checkpts,
                   const std::vector<int>& checkpts_sero );

/*
 Simulate the spread of DENV (no ZIKV). Yields the age distribution of a single Dengue serotype and of all serotypes together between a set of user-defined checkpoints, distinguishing by infection parity
 */
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
                   const std::vector<int>& checkpts );



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
                   const std::vector<int>& checkpts_sero );

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
                       const double& sero_a0 );

/*
 Simulate the spread of DENV (no ZIKV). Remove a single DENV strain for some time and then
 introduce it back.
*/
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
                   const std::vector<int>& checkpts );



/*
 Simulate the dynamics of ZIKV
 */
std::vector<std::vector<std::vector<std::vector<int>>>> simulate_ZIKV(const int& Nsim,
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
                                                                      const int& nZIKV0 = 10 );

/*
 Simulate the dynamics of ZIKV
 Uses realistic seroprevalence profile
 Measures incidence by age
 */
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
                                                                      const int& nZIKV0 = 10 );


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
                                                                      const int& nZIKV0 = 10 );
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
                                                  const double& chiZD2 );






#endif /* pyFlavi_h */
