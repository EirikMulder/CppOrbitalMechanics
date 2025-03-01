#include <iostream>
// #include <array>
#include <memory>
#include <functional>
#include <vector>
#include <matplot/matplot.h>
#include <cmath>
#include <algorithm>

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

    std::string display() const {
      /*return "<" + x + ", " + y + ", " + z + ", " + vx + ", " + vy + ", " + vz + ">";*/
      return std::format("<{}, {}, {}, {}, {}, {}>", x, y, z, vx, vy, vz);
    }
};

StateVector euler_step(std::function<StateDerivative(double, StateVector)> f, double t, const StateVector &y, double dt) {
  return y + dt * f(dt, y);
}

std::unique_ptr<std::vector<StateVector>> euler_propagation(std::function<StateDerivative(double, StateVector)> f, double t0, const StateVector &y0, double dt, double tf) {
  StateVector y = y0;
  std::unique_ptr<std::vector<StateVector>> states(new std::vector<StateVector>());
  for (double t = t0; t < tf; t += dt) {
    states->push_back(y);
    y = euler_step(f, t, y, dt);
  }
  states->push_back(y);
  return states;
}

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

void plot_orbit(std::unique_ptr<std::vector<StateVector>> const &states) {
  using namespace matplot;
  std::vector<StateVector> st = *states;
  std::vector<double> xs(states->size());
  std::vector<double> ys(states->size());
  std::vector<double> zs(states->size());
  std::transform(states->begin(), states->end(), xs.begin(), [] (const StateVector &i) {return i.x;});
  std::transform(states->begin(), states->end(), ys.begin(), [] (const StateVector &i) {return i.y;});
  std::transform(states->begin(), states->end(), zs.begin(), [] (const StateVector &i) {return i.z;});

  std::cout << "states size = " << states->size() << std::endl;
  std::cout << "xs size = " << xs.size() << std::endl;
  std::cout << "ys size = " << ys.size() << std::endl;
  std::cout << "zs size = " << zs.size() << std::endl;

  auto p = plot3(xs, ys, zs);
  p->parent()->x_axis().tick_label_format("s");
  xlabel("X");
  ylabel("Y");
  zlabel("Z");
  grid(on);
  show();
  /*std::cout << "Plotting!" << std::endl;*/
}

void save_orbit(std::unique_ptr<std::vector<StateVector>> const &states, std::string filename) {
  
}


int main() {
  std::cout << "Hello World!" << std::endl;
  double altitude = 400;
  double radius = R_E + altitude;
  double vel = sqrt(MU / radius);
  // std::array<double, 6> y0 = {radius, 0, 0, 0, vel, 0};
  StateVector y0 = StateVector(radius, 0, 0, 0, vel, 0);
  std::cout << "Initial Radius: " << radius << " [km], Initial Velocity: "
            << vel << " [km/s]" << std::endl;
  double dt = 1;
  double tf = 10000;
  auto states = euler_propagation(keplerian_dynamics, 0, y0, dt, tf);
  std::cout << "Propagation Finished!" << std::endl;
  /*for (StateVector y : *states) {*/
  /*  std::cout << y.display() << "\n";*/
  /*}*/
  /*std::cout << std::endl;*/
  plot_orbit(states);
}
