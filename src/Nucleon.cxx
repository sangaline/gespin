#include "Nucleon.h"

#include "math.h"

Nucleon::Nucleon(double x, double y, double z) : x(x), y(y), z(z) {

}

double Nucleon::R() {
    return sqrt(x*x + y*y + z*z);
}

double Nucleon::Phi() {
    return atan2(y, x);
}

double Nucleon::Theta() {
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
