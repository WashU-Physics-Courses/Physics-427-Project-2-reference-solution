#include "fluid_solver_2d.h"

int main() {
  int Nx = 258;
  int Ny = 258;
  double Lx = 1.0;
  double Ly = 1.0;
  double gamma = 5.0 / 3.0;
  fluid_solver_2d solver(Lx, Ly, Nx, Ny, gamma, 0.01);

  // initial conditions
  std::vector<double> rho0(Nx * Ny);
  std::vector<double> vx0(Nx * Ny);
  std::vector<double> vy0(Nx * Ny);
  std::vector<double> P0(Nx * Ny);

  double sigma = 0.05;
  for (int j = 0; j < Ny; j++) {
    for (int i = 0; i < Nx; i++) {
      int idx = i + j * Nx;
      double x = (i - 1) * Lx / (Nx - 2);
      double y = (j - 1) * Ly / (Ny - 2);

      // if (y < 0.7 && y > 0.3 && x < 0.7 && x > 0.3) {
      if (y < 0.55 && y > 0.45) {
        rho0[idx] = 4.0;
        vx0[idx] = 0.5;
      } else {
        rho0[idx] = 1.0;
        vx0[idx] = -0.5;
      }
      vy0[idx] = 0.01 * (std::sin(6.0 * M_PI * x)) *
                 (std::exp(-std::pow((y - 0.75) / sigma, 2.0)) +
                  std::exp(-std::pow((y - 0.25) / sigma, 2.0)));
      P0[idx] = 2.5;
    }
  }

  solver.init(rho0, vx0, vy0, P0);
  solver.solve(0.0, 4.0);

  return 0;
}
