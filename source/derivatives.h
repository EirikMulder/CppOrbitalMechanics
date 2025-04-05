#pragma once
#include <functional>
#include <memory>
#include <vector>

namespace orbital {

  namespace coefficients {
inline constexpr std::vector<std::vector<double>> b_table = {
    {0, 0, 0, 0, 0},
    {2. / 9, 0, 0, 0, 0},
    {1. / 12, 1. / 4, 0, 0, 0},
    {69. / 128, -243. / 128, 135. / 64, 0, 0},
    {-17. / 12, 27. / 4, -27. / 5, 16. / 15, 0},
    {65. / 432, -5. / 16, 13. / 16, 4. / 27, 5. / 144},
};

inline constexpr std::vector<double> a_table = {0, 2. / 9, 1. / 3, 3. / 4, 1, 5. / 6};

inline constexpr std::vector<double> c_table = {1. / 9, 0, 9. / 20, 16. / 45, 1. / 12};

inline constexpr std::vector<double> ch_table = {47. / 450, 0,       12. / 25,
                                      32. / 225, 1. / 30, 6. / 25};

inline constexpr std::vector<double> ct_table = {1. / 150, 0,       -3. / 100,
                                      16. / 75, 1. / 20, -6. / 25};

  }

class StateDerivative {
public:
  double vx;
  double vy;
  double vz;
  double ax;
  double ay;
  double az;

  StateDerivative(double vx, double vy, double vz, double ax, double ay,
                  double az);
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

struct DynamicStepResult {
  StateVector y;
  double h_new;
};

double abs(const StateDerivative &derivative);

template <typename T>
StateVector euler_step(std::function<StateDerivative(double, StateVector, T)> f,
                       double t, const StateVector &y, double dt, T cfg) {
  return y + dt * f(t, y, cfg);
}

template <typename T>
StateVector rk4_step(std::function<StateDerivative(double, StateVector, T)> f,
                     double t, const StateVector &y, double dt, T cfg) {
  auto k1 = dt * f(t, y, cfg);
  auto k2 = dt * f(t + dt / 2, y + 1. / 2 * k1, cfg);
  auto k3 = dt * f(t + dt / 2, y + 1. / 2 * k2, cfg);
  auto k4 = dt * f(t + dt, y + k3, cfg);
  return y + (1. / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
}

template <typename T>
DynamicStepResult
rk45_step(std::function<StateDerivative(double, StateVector, T)> f, double t,
          const StateVector &y, double h, double epsilon, T cfg) {
  double truncation_error;
  StateVector y_new = y;
  do {
    auto k1 = h * f(t + a_table[0] * h, y, cfg);
    auto k2 = h * f(t + a_table[1] * h, y + b_table[1][0] * k1, cfg);
    auto k3 = h * f(t + a_table[2] * h,
                    y + b_table[2][0] * k1 + b_table[2][1] * k2, cfg);
    auto k4 =
        h * f(t + a_table[3] * h,
              y + b_table[3][0] * k1 + b_table[3][1] * k2 + b_table[3][2] * k3,
              cfg);
    auto k5 = h * f(t + a_table[4] * h,
                    y + b_table[4][0] * k1 + b_table[4][1] * k2 +
                        b_table[4][2] * k3 + b_table[4][3] * k4,
                    cfg);
    auto k6 =
        h * f(t + a_table[5] * h,
              y + b_table[5][0] * k1 + b_table[5][1] * k2 + b_table[5][2] * k3 +
                  b_table[5][3] * k4 + b_table[5][4] * k5,
              cfg);
    y_new = y + ch_table[0] * k1 + ch_table[1] * k2 + ch_table[2] * k3 +
            ch_table[3] * k4 + ch_table[4] * k5 + ch_table[5] * k6;
    truncation_error =
        abs(ct_table[0] * k1 + ct_table[1] * k2 + ct_table[2] * k3 +
            ct_table[3] * k4 + ct_table[4] * k5 + ct_table[5] * k6);
    h = 0.9 * h * pow(epsilon / truncation_error, 1. / 5);
  } while (truncation_error > epsilon);
  return DynamicStepResult{y_new, h};
}

template <typename T>
std::vector<StateVector>
euler_propagation(std::function<StateDerivative(double, StateVector, T)> f,
                  double t0, const StateVector &y0, double dt, double tf,
                  T cfg) {
  StateVector y = y0;
  // std::unique_ptr<std::vector<StateVector>> states(
  //     new std::vector<StateVector>());
  std::vector<StateVector> states;
  int expected_steps = ceil(tf / dt);
  states.reserve(expected_steps);
  for (double t = t0; t < tf; t += dt) {
    states.push_back(y);
    y = euler_step(f, t, y, dt, cfg);
  }
  states.push_back(y);
  return states;
}

template <typename T>
std::vector<StateVector>
rk4_propagation(std::function<StateDerivative(double, StateVector)> f,
                double t0, const StateVector &y0, double dt, double tf, T cfg) {
  StateVector y = y0;
  std::vector<StateVector> states;
  int expected_steps = ceil(tf / dt);
  states.reserve(expected_steps);
  for (double t = t0; t < tf; t += dt) {
    states.push_back(y);
    y = rk4_step(f, t, y, dt, cfg);
  }
  states.push_back(y);
  return states;
}

template <typename T>
std::vector<StateVector>
rk45_propagation(std::function<StateDerivative(double, StateVector, T)> f,
                 double t0, const StateVector &y0, double dt, double tf,
                 double epsilon, T cfg) {
  StateVector y = y0;
  std::vector<StateVector> states;
  int expected_steps = ceil(tf / dt);
  states.reserve(expected_steps);
  double t = t0;
  while (t < tf) {
    DynamicStepResult result = rk45_step(f, t, y, dt, epsilon, cfg);
    y = result.y;
    dt = result.h_new;
    states.push_back(y);
    t += dt;
  }
  states.push_back(y);
  return states;
}

bool save_orbit(std::vector<StateVector> const &states,
                std::string filename);

}
