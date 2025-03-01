#include <iostream>
#include <fstream>
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

    StateDerivative operator+(StateDerivative const &b) const {
      return StateDerivative(vx + b.vx, vy + b.vy, vz + b.vz, ax + b.ax, ay + b.ay, az + b.az);
    }
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
  return y + dt * f(t, y);
}

StateVector rk4_step(std::function<StateDerivative(double, StateVector)> f, double t, const StateVector &y, double dt) {
  auto k1 = dt * f(t, y);
  auto k2 = dt * f(t + dt/2, y + 1./2 * k1);
  auto k3 = dt * f(t + dt/2, y + 1./2 * k2);
  auto k4 = dt * f(t + dt, y + k3);
  return y + (1./6) * (k1 + 2*k2 + 2*k3 + k4);
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
}
