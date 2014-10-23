#pragma once

#include <vector>
#include <functional>

class Nucleon;

class NucleonCollection {
  public:
    typedef std::vector<Nucleon*> nucleon_array;
  private:
    nucleon_array ***cubes;
    nucleon_array ordered;

    double likelihood;
    unsigned int nucleon_count;
    double LikelihoodContribution(Nucleon *nucleon);

    //these specify one dimension of one octet
    //the full number of cubes is thus (2*units)^3
    unsigned int units;
    double length, cube_length;
    int pairwise_units;

    void RemoveNucleonFromCube(Nucleon *nucleon);
    void RemoveNucleonFromOrdered(Nucleon *nucleon);

    std::function<double(Nucleon&)> single_likelihood;
    std::function<double(Nucleon&, Nucleon&)> pairwise_likelihood;
  public:
    NucleonCollection(double pairwise_max = 0, unsigned int units = 10, double length = 20);
    NucleonCollection(const NucleonCollection &nucleon_collection);
    ~NucleonCollection();

    unsigned int AddNucleon(const Nucleon &nucleon);
    void Reset(bool delete_nucleons = true);

    Nucleon& operator [](int i) { return *ordered[i]; }
    const Nucleon& operator [](int i) const { return *ordered[i]; }
    void RemoveNucleon(Nucleon *nucleon);
    unsigned int NucleonCount() const { return nucleon_count; }
    unsigned int Units() { return units; }
    double Length() { return length; }

    void SetNucleonX(Nucleon *nucleon, double x);
    void SetNucleonY(Nucleon *nucleon, double y);
    void SetNucleonZ(Nucleon *nucleon, double z);
    void SetNucleonPosition(Nucleon *nucleon, double x, double y, double z);

    nucleon_array& FindCube(double x, double y, double z, int &i, int &j, int &k) const;
    nucleon_array& FindCube(double x, double y, double z) const;
    nucleon_array& FindCube(const Nucleon &nucleon) const;

    virtual double SingleLikelihood(Nucleon &nucleon) { return single_likelihood(nucleon); }
    virtual double PairwiseLikelihood(Nucleon &nucleon1, Nucleon &nucleon2) { return pairwise_likelihood(nucleon1, nucleon2); }
    double Likelihood() { return likelihood; }
};
