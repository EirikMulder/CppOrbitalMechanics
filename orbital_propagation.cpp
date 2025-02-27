#include <iostream>
#include <array>
#include <memory>

const double R_E = 6378;
const double MU = 398600;


int main() {
  std::cout << "Hello World!" << std::endl;
  double altitude = 400;
  double radius = R_E + altitude;
  double vel = sqrt(MU / radius);
  std::array<double, 6> y0 = {radius, 0, 0, 0, vel, 0};
  std::cout << "Initial Radius: " << radius << " [km], Initial Velocity: "
            << vel << " [km/s]" << std::endl;
}
