#pragma once

class Nucleon {
  private:
    double x, y, z;
  public:
    Nucleon(double x = 0, double y = 0, double z = 0);

    double X() { return x; }
    double Y() { return y; }
    double Z() { return z; }
    void SetX(double newx) { x = newx; }
    void SetY(double newy) { y = newy; }
    void SetZ(double newz) { z = newz; }
    void SetPosition(double newx, double newy, double newz) { x = newx; y = newy; z = newz; }
};
