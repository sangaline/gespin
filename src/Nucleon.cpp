#include "Nucleon.h"
#include "NucleonCollection.h"

#include "math.h"

Nucleon::Nucleon(double x, double y, double z, double radius, NucleonIdentity identity)
  : x(x), y(y), z(z), radius(radius), identity(identity), parent(0) {

}

//the parent is maintained and a new cube is calculated
Nucleon& Nucleon::operator=(const Nucleon &nucleon) {
    x = nucleon.x;
    y = nucleon.y;
    z = nucleon.z;
    radius = nucleon.radius;
    identity = nucleon.identity;
    cube = &parent->FindCube(x, y, z);

    return *this;
}

double Nucleon::R()  const {
    return sqrt(x*x + y*y + z*z);
}

double Nucleon::Phi() const {
    return atan2(y, x);
}

double Nucleon::Theta() const {
    return atan2(sqrt(x*x + y*y), z);
}

void Nucleon::SetR(double new_r) {
    const double factor = new_r/R();
    SetPosition(factor*X(), factor*Y(), factor*Z());
}

void Nucleon::SetTheta(double new_theta) {
    const double r = R();
    const double phi = Phi();
    const double sin_theta = sin(new_theta);
    SetX(r*cos(phi)*sin_theta);
    SetY(r*sin(phi)*sin_theta);
    SetZ(r*cos(new_theta));
}

void Nucleon::SetPhi(double new_phi) {
    const double r_transverse = sqrt(X()*X() + Y()*Y());
    SetX(r_transverse*cos(new_phi));
    SetY(r_transverse*sin(new_phi));
}
void Nucleon::SetPosition(double new_x, double new_y, double new_z) { 
    if(parent) {
        parent->SetNucleonPosition(this, new_x, new_y, new_z);
    }
    else {
        x = new_x;
        y = new_y;
        z = new_z;
    }
}
