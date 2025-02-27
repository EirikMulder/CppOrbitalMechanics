#include <iostream>
#include <array>
#include <memory>

const double R_E = 6378;
const double MU = 398600;

class StateDerivative {
  public:
    double vx;
    double vy;
    double vz;
    double ax;
    double ay;
    double az;

    StateDerivative(double vx, double vy, double vz, double ax, double ay, double az)
    : vx(vx), vy(vy), vz(vz), ax(ax), ay(ay), az(az) { }
};

StateDerivative operator*(double a, StateDerivative const &b) {
  return StateDerivative(b.vx * a, b.vy * a, b.vz * a, b.ax * a, b.ay * a, b.az * a);
}

class StateVector {
  public:
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;

    StateVector(double x, double y, double z, double vx, double vy, double vz)
    : x(x), y(y), z(z), vx(vx), vy(vy), vz(vz) { }

    StateVector operator+(StateDerivative const &b) const {
      return StateVector(x + b.vx, y + b.vy, z + b.vz, vx + b.ax, vy + b.ay, vz + b.az);
    }
};


int main() {
  std::cout << "Hello World!" << std::endl;
  double altitude = 400;
  double radius = R_E + altitude;
  double vel = sqrt(MU / radius);
  // std::array<double, 6> y0 = {radius, 0, 0, 0, vel, 0};
  StateVector y0 = StateVector(radius, 0, 0, 0, vel, 0);
  std::cout << "Initial Radius: " << radius << " [km], Initial Velocity: "
            << vel << " [km/s]" << std::endl;
}
