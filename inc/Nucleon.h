#pragma once

enum class NucleonIdentity { unspecified, proton, neutron, antineutron, antiproton };
class Nucleon {
  private:
    double x, y, z;
    double radius;
    NucleonIdentity identity;
  public:
    Nucleon(double x = 0, double y = 0, double z = 0, double radius = 0.8768, NucleonIdentity identity = NucleonIdentity::unspecified);

    //property getters
    double Radius() { return radius; }
    NucleonIdentity Identity() { return identity; }

    //property setters
    void SetRadius(double new_radius) { radius = new_radius; }
    void SetIdentity(NucleonIdentity new_identity) { identity = new_identity; }

    //position getters
    double X() { return x; }
    double Y() { return y; }
    double Z() { return z; }
    double R();
    double Theta();
    double Phi();

    //position setters
    void SetX(double new_x) { x = new_x; }
    void SetY(double new_y) { y = new_y; }
    void SetZ(double new_z) { z = new_z; }
    void SetR(double new_r);
    void SetTheta(double new_theta);
    void SetPhi(double new_phi);
    void SetPosition(double new_x, double new_y, double new_z) { x = new_x; y = new_y; z = new_z; }
};
