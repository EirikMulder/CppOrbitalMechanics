#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "derivatives.h"

// Define absolute value of a vector
double abs(const StateDerivative &derivative) {
  return sqrt(pow(derivative.vx, 2) + pow(derivative.vy, 2) +
              pow(derivative.vz, 2) + pow(derivative.ax, 2) +
              pow(derivative.ay, 2) + pow(derivative.az, 2));
}

StateDerivative::StateDerivative(double vx, double vy, double vz, double ax,
                                 double ay, double az)
    : vx(vx), vy(vy), vz(vz), ax(ax), ay(ay), az(az) {}

StateDerivative StateDerivative::operator+(StateDerivative const &b) const {
  return StateDerivative(vx + b.vx, vy + b.vy, vz + b.vz, ax + b.ax, ay + b.ay,
                         az + b.az);
}

StateDerivative operator*(double a, StateDerivative const &b) {
  return StateDerivative(b.vx * a, b.vy * a, b.vz * a, b.ax * a, b.ay * a,
                         b.az * a);
}

StateVector::StateVector(double x, double y, double z, double vx, double vy,
                         double vz)
    : x(x), y(y), z(z), vx(vx), vy(vy), vz(vz) {}

StateVector StateVector::operator+(StateDerivative const &b) const {
  return StateVector(x + b.vx, y + b.vy, z + b.vz, vx + b.ax, vy + b.ay,
                     vz + b.az);
}

std::string StateVector::display() const {
  /*return "<" + x + ", " + y + ", " + z + ", " + vx + ", " + vy + ", " + vz +
   * ">";*/
  return std::format("<{}, {}, {}, {}, {}, {}>", x, y, z, vx, vy, vz);
}

void save_orbit(std::unique_ptr<std::vector<StateVector>> const &states,
                std::string filename) {
  std::ofstream file;
  file.open(filename);
  file << "x,y,z,vx,vy,vz";
  for (StateVector i : *states) {
    file << "\n";
    file << i.x << "," << i.y << "," << i.z << ",";
    file << i.vx << "," << i.vy << "," << i.vz;
  }
  file << std::endl;
}
