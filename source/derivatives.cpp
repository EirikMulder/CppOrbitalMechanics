#include <iostream>
#include <fstream>
#include <memory>
#include <functional>
#include <vector>
#include <cmath>

#include "derivatives.h"

const std::vector<std::vector<double>> b_table = {
  {0, 0, 0, 0, 0},
  {2./9, 0, 0, 0, 0},
  {1./12, 1./4, 0, 0, 0},
  {69./128, -243./128, 135./64, 0, 0},
  {-17./12, 27./4, -27./5, 16./15, 0},
  {65./432, -5./16, 13./16, 4./27, 5./144},
};

const std::vector<double> a_table = {0, 2./9, 1./3, 3./4, 1, 5./6};

const std::vector<double> c_table = {1./9, 0, 9./20, 16./45, 1./12};

const std::vector<double> ch_table = {47./450, 0, 12./25, 32./225, 1./30, 6./25};

const std::vector<double> ct_table = {1./150, 0, -3./100, 16./75, 1./20, -6./25};


StateDerivative::StateDerivative(double vx, double vy, double vz, double ax, double ay, double az)
    : vx(vx), vy(vy), vz(vz), ax(ax), ay(ay), az(az) { }

StateDerivative StateDerivative::operator+(StateDerivative const &b) const {
      return StateDerivative(vx + b.vx, vy + b.vy, vz + b.vz, ax + b.ax, ay + b.ay, az + b.az);
    }

StateDerivative operator*(double a, StateDerivative const &b) {
  return StateDerivative(b.vx * a, b.vy * a, b.vz * a, b.ax * a, b.ay * a, b.az * a);
}

StateVector::StateVector(double x, double y, double z, double vx, double vy, double vz)
: x(x), y(y), z(z), vx(vx), vy(vy), vz(vz) { }

StateVector StateVector::operator+(StateDerivative const &b) const {
  return StateVector(x + b.vx, y + b.vy, z + b.vz, vx + b.ax, vy + b.ay, vz + b.az);
}

std::string StateVector::display() const {
  /*return "<" + x + ", " + y + ", " + z + ", " + vx + ", " + vy + ", " + vz + ">";*/
  return std::format("<{}, {}, {}, {}, {}, {}>", x, y, z, vx, vy, vz);
}

struct DynamicStepResult {
  StateVector y;
  double h_new;
};

StateVector euler_step(std::function<StateDerivative(double, StateVector)> f, double t, const StateVector &y, double dt) {
  return y + dt * f(t, y);
}

StateVector rk4_step(std::function<StateDerivative(double, StateVector)> f, double t, const StateVector &y, double dt) {
  auto k1 = dt * f(t, y);
  auto k2 = dt * f(t + dt/2, y + 1./2 * k1);
  auto k3 = dt * f(t + dt/2, y + 1./2 * k2);
  auto k4 = dt * f(t + dt, y + k3);
  return y + (1./6) * (k1 + 2*k2 + 2*k3 + k4);
}

StateVector rk45_step(std::function<StateDerivative(double, StateVector)> f, double t, const StateVector &y, double h) {

}

std::unique_ptr<std::vector<StateVector>> euler_propagation(std::function<StateDerivative(double, StateVector)> f, double t0, const StateVector &y0, double dt, double tf) {
  StateVector y = y0;
  std::unique_ptr<std::vector<StateVector>> states(new std::vector<StateVector>());
  int expected_steps = ceil(tf / dt);
  states->reserve(expected_steps);
  for (double t = t0; t < tf; t += dt) {
    states->push_back(y);
    y = euler_step(f, t, y, dt);
  }
  states->push_back(y);
  return states;
}

std::unique_ptr<std::vector<StateVector>> rk4_propagation(std::function<StateDerivative(double, StateVector)> f, double t0, const StateVector &y0, double dt, double tf) {
  StateVector y = y0;
  std::unique_ptr<std::vector<StateVector>> states(new std::vector<StateVector>());
  int expected_steps = ceil(tf / dt);
  states->reserve(expected_steps);
  for (double t = t0; t < tf; t += dt) {
    states->push_back(y);
    y = rk4_step(f, t, y, dt);
  }
  states->push_back(y);
  return states;
}


void save_orbit(std::unique_ptr<std::vector<StateVector>> const &states, std::string filename) {
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
