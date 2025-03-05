#include <iostream>
#include <functional>
#include <cmath>

#include "derivatives.h"

const double R_E = 6378;
const double MU = 398600;

StateDerivative keplerian_dynamics(double t, StateVector const &y) {
  double radius = sqrt(pow(y.x, 2) + pow(y.y, 2) + pow(y.z, 2));
  double acc_coefficient = -MU / pow(radius, 3);
  return StateDerivative(
    y.vx,
    y.vy,
    y.vz,
    acc_coefficient * y.x,
    acc_coefficient * y.y,
    acc_coefficient * y.z
  );
}



int main(int argc, char *argv[]) {
  double altitude = 400;
  double radius = R_E + altitude;
  double vel = sqrt(MU / radius);
  double x0 = radius;
  double vy0 = vel;
  double y0 = 0, z0 = 0;
  double vx0 = 0, vz0 = 0;
  // CLI Arguments
  std::string filename = "orbit.csv";
  switch (argc) {
    case 1:
      break;
    case 2:
      filename = argv[1];
      break;
    case 8:
      filename = argv[7];
    case 7:
      x0 = std::stod(argv[1]);
      y0 = std::stod(argv[2]);
      z0 = std::stod(argv[3]);
      vx0 = std::stod(argv[4]);
      vy0 = std::stod(argv[5]);
      vz0 = std::stod(argv[6]);
      std::cout << std::format("r = < {}, {}, {} >, v = < {}, {}, {} >", x0, y0, z0, vx0, vy0, vz0) << std::endl;
      break;
    default:
      std::cout << "Invalid number of arguments!\n";
      std::cout << "\tOptions:\n";
      std::cout << "\t<filename>\n";
      std::cout << "\t<x0> <y0> <z0> <vx0> <vy0> <vz0>\n";
      std::cout << "\t<x0> <y0> <z0> <vx0> <vy0> <vz0> <filename>" << std::endl;
  }

  std::cout << "Hello World!" << std::endl;
  // std::array<double, 6> y0 = {radius, 0, 0, 0, vel, 0};
  StateVector Y0 = StateVector(x0, y0, z0, vx0, vy0, vz0);
  std::cout << "Initial Radius: " << radius << " [km], Initial Velocity: "
            << vel << " [km/s]" << std::endl;
  double dt = 1;
  double tf = 10000;
  auto states = rk4_propagation(keplerian_dynamics, 0, Y0, dt, tf);
  std::cout << "Propagation Finished!" << std::endl;
  save_orbit(states, filename);
  std::cout << "Saved Orbital Data to " << filename << std::endl;
  /*plot_orbit(states);*/
}
