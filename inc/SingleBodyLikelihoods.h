#pragma once

#include <functional>

class Nucleon;
namespace SingleBodyLikelihoods {
    class WoodsSaxon {
      public:
          double R, a, B_20, B_40, w;
          WoodsSaxon(double R, double a, double B_20 = 0, double B_40 = 0, double w = 0);
          bool Deconvolve(double sigma);

          double Likelihood(const double r, const double theta) const;
          double Likelihood(const Nucleon& nucleon) const;
          static double Likelihood(const double r, const double theta, const double R, const double a, const double B_20 = 0, const double B_40 = 0, const double w = 0);
          std::function<double(const Nucleon&)> LikelihoodFunction() const;

    };
    typedef std::function<double(Nucleon& nucleon)> single_body_function;

    WoodsSaxon Au197();
    WoodsSaxon Pb208();
    WoodsSaxon Cu63();
    WoodsSaxon U238();
}
