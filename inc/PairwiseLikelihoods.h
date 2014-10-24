#pragma once

#include <functional>

class Nucleon;
namespace PairwiseLikelihoods {
    class GaussianRepulsion {
        double beta, strength;
      public:
          GaussianRepulsion(double sigma = 0.74536, double strength = 1);
          double Likelihood(double r) const;
          double Likelihood(const Nucleon& n1, const Nucleon& n2) const;
          static double Likelihood(const Nucleon& n1, const Nucleon& n2, double beta, double strength = 1);
          std::function<double(const Nucleon&, const Nucleon&)> LikelihoodFunction() const;

          double Beta() { return beta; }
          double Sigma();
          double Strength() { return strength; }

          bool SetBeta(double new_beta);
          bool SetSigma(double new_sigma);
          bool SetStrength(double new_strength);
    };
};
