#include "PairwiseLikelihoods.h"
#include "Nucleon.h"

#include "math.h"

PairwiseLikelihoods::GaussianRepulsion::GaussianRepulsion(double sigma, double strength) {
    SetSigma(sigma);
    SetStrength(strength);
}

double PairwiseLikelihoods::GaussianRepulsion::Likelihood(double r) const {
    return 1 - strength*exp(-beta*r*r);
}

double PairwiseLikelihoods::GaussianRepulsion::Likelihood(const Nucleon& nucleon1, const Nucleon& nucleon2) const {
    return Likelihood(nucleon1, nucleon2, beta, strength);
}

double PairwiseLikelihoods::GaussianRepulsion::Likelihood(const Nucleon& nucleon1, const Nucleon& nucleon2, double beta, double strength) {
    const double r2 = pow(nucleon1.X() - nucleon2.X(), 2)
                    + pow(nucleon1.Y() - nucleon2.Y(), 2)
                    + pow(nucleon1.Z() - nucleon2.Z(), 2);
    return 1 - strength*exp(-beta*r2);
}

std::function<double(const Nucleon&, const Nucleon&)> PairwiseLikelihoods::GaussianRepulsion::LikelihoodFunction() const {
    const double beta_ = beta, strength_ = strength;
    return [beta_, strength_](const Nucleon& nucleon1, const Nucleon& nucleon2) {
        return Likelihood(nucleon1, nucleon2, beta_, strength_);
    };
}

double PairwiseLikelihoods::GaussianRepulsion::Sigma() { 
    return sqrt(1/(2*beta));
}

bool PairwiseLikelihoods::GaussianRepulsion::SetBeta(double new_beta) {
    if(new_beta < 0) {
        return false;
    }
    beta = new_beta;
    return true;
}

bool PairwiseLikelihoods::GaussianRepulsion::SetSigma(double new_sigma) {
    if(new_sigma < 0) {
        return false;
    }
    beta = 0.5/(new_sigma*new_sigma);
    return true;
}

bool PairwiseLikelihoods::GaussianRepulsion::SetStrength(double new_strength) {
    strength = new_strength;
    if(strength > 1) {
        strength = 1;
        return false;
    }
    return true;
}
