#include "SingleBodyLikelihoods.h"
#include <Nucleon.h>

#include <boost/math/special_functions/spherical_harmonic.hpp>
#include "math.h"

using namespace boost::math;

SingleBodyLikelihoods::WoodsSaxon::WoodsSaxon(double R, double a, double B_20, double B_40, double w) 
  : R(R), a(a), B_20(B_20), B_40(B_40), w(w) { }

bool SingleBodyLikelihoods::WoodsSaxon::Deconvolve(double sigma) {
    const double sigma_term = sigma*sigma*M_PI/8;
    const double a2 = a*a;
    //deconvolution is not possible
    if(a2 > sigma_term) {
        return false;
    }
    a = sqrt(sigma_term - a2);
    return true;
}

double SingleBodyLikelihoods::WoodsSaxon::Likelihood(const double r, const double theta) const {
    return Likelihood(r, theta, R, a, B_20, B_40, w);
}

double SingleBodyLikelihoods::WoodsSaxon::Likelihood(const Nucleon& nucleon) const {
    return Likelihood(nucleon.R(), nucleon.Theta());
}

double SingleBodyLikelihoods::WoodsSaxon::Likelihood(const double r, const double theta, const double R, const double a, const double B_20, const double B_40, const double w) {
    double deformation = 1;
    if(B_20 != 0) {
        deformation += B_20*spherical_harmonic_r(2, 0, theta, 0);
    }
    if(B_40 != 0) {
        deformation += B_40*spherical_harmonic_r(4, 0, theta, 0);
    }
    const double x = (r - R*deformation)/(a*deformation);
    double w_effective = w;
    if(w_effective != 0 && r < R) {
        w_effective *= pow(r/R, 2);
    }

    return (1 + w_effective)/(1 + exp(x));
}

std::function<double(const Nucleon&)> SingleBodyLikelihoods::WoodsSaxon::LikelihoodFunction() const {
    const double R_ = R, a_ = a, B_20_ = B_20, B_40_ = B_20, w_ = w;
    return [R_, a_, B_20_, B_40_, w_](const Nucleon &nucleon) {
        return SingleBodyLikelihoods::WoodsSaxon::Likelihood(nucleon.R(), nucleon.Theta(), R_, a_, B_20_, B_40_, w_);
    };
}

SingleBodyLikelihoods::WoodsSaxon SingleBodyLikelihoods::Au197() {
    return SingleBodyLikelihoods::WoodsSaxon(6.38, 0.535, -0.131, -0.031);
}

SingleBodyLikelihoods::WoodsSaxon SingleBodyLikelihoods::Pb208() {
    return SingleBodyLikelihoods::WoodsSaxon(6.624, 0.549);
}

SingleBodyLikelihoods::WoodsSaxon SingleBodyLikelihoods::Cu63() {
    return SingleBodyLikelihoods::WoodsSaxon(4.24, 0.586, 0.162, -0.006);
}

SingleBodyLikelihoods::WoodsSaxon SingleBodyLikelihoods::U238() {
    return SingleBodyLikelihoods::WoodsSaxon(6.8, 0.5, 0.254, 0.052);
}
