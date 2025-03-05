#pragma once

#include <memory>
#include <functional>
#include <vector>

class StateDerivative {
  public:
    double vx;
    double vy;
    double vz;
    double ax;
    double ay;
    double az;

    StateDerivative(double vx, double vy, double vz, double ax, double ay, double az);
    StateDerivative operator+(StateDerivative const &b) const;
};

StateDerivative operator*(double a, StateDerivative const &b);

class StateVector {
  public:
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;

    StateVector(double x, double y, double z, double vx, double vy, double vz);
    StateVector operator+(StateDerivative const &b) const;
    std::string display() const;
};

StateVector euler_step(std::function<StateDerivative(double, StateVector)> f, double t, const StateVector &y, double dt);

StateVector rk4_step(std::function<StateDerivative(double, StateVector)> f, double t, const StateVector &y, double dt);

std::unique_ptr<std::vector<StateVector>> euler_propagation(std::function<StateDerivative(double, StateVector)> f, double t0, const StateVector &y0, double dt, double tf);

std::unique_ptr<std::vector<StateVector>> rk4_propagation(std::function<StateDerivative(double, StateVector)> f, double t0, const StateVector &y0, double dt, double tf);

void save_orbit(std::unique_ptr<std::vector<StateVector>> const &states, std::string filename);

