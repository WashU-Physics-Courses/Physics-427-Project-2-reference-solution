#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

double square(double x) { return x * x; }

class fluid_solver_2d {
public:
  fluid_solver_2d(double Lx0, double Ly0, int Nx0, int Ny0, double gamma0,
                  double output_dt0)
      : gamma(gamma0), Lx(Lx0), Ly(Ly0), Nx(Nx0), Ny(Ny0) {
    output_dt = output_dt0;
    dx = Lx / (Nx - 2);
    dy = Ly / (Ny - 2);

    rho.resize(Nx * Ny);
    vx.resize(Nx * Ny);
    vy.resize(Nx * Ny);
    P.resize(Nx * Ny);

    mass.resize(Nx * Ny);
    mom_x.resize(Nx * Ny);
    mom_y.resize(Nx * Ny);
    energy.resize(Nx * Ny);

    rho_tmp.resize(Nx * Ny);
    vx_tmp.resize(Nx * Ny);
    vy_tmp.resize(Nx * Ny);
    P_tmp.resize(Nx * Ny);

    rho_Lx.resize(Nx * Ny);
    rho_Rx.resize(Nx * Ny);
    rho_Ly.resize(Nx * Ny);
    rho_Ry.resize(Nx * Ny);

    vx_Lx.resize(Nx * Ny);
    vx_Rx.resize(Nx * Ny);
    vx_Ly.resize(Nx * Ny);
    vx_Ry.resize(Nx * Ny);
    vy_Lx.resize(Nx * Ny);
    vy_Rx.resize(Nx * Ny);
    vy_Ly.resize(Nx * Ny);
    vy_Ry.resize(Nx * Ny);
    P_Lx.resize(Nx * Ny);
    P_Rx.resize(Nx * Ny);
    P_Ly.resize(Nx * Ny);
    P_Ry.resize(Nx * Ny);

    mass_flux_x.resize(Nx * Ny);
    mass_flux_y.resize(Nx * Ny);
    momx_flux_x.resize(Nx * Ny);
    momx_flux_y.resize(Nx * Ny);
    momy_flux_x.resize(Nx * Ny);
    momy_flux_y.resize(Nx * Ny);
    energy_flux_x.resize(Nx * Ny);
    energy_flux_y.resize(Nx * Ny);
  }

  ~fluid_solver_2d() {}

  void primitive_to_conserved() {
    // TODO: Compute conserved variables from primitive ones
  }

  void conserved_to_primitive() {
    // TODO: Compute primitive variables from conserved ones
  }

  void init(const std::vector<double> &rho0, const std::vector<double> &vx0,
            const std::vector<double> &vy0, const std::vector<double> &P0) {
    // TODO: Initialize the primitive variables using the given initial
    // condition
  }

  double find_dt() {
    // TODO: Find the optimal dt that satisfies the CFL condition, and return
    // its value
  }

  void solve(double t0, double t_end) {
    // Solve the fluid equations, starting from t0 and stoping at t_end
    double t = t0;
    int n = 0; // n labels the output file
    while (t < t_end) {
      if (t >= output_dt * n) {
        output(n);
        n += 1;
      }
      double dt = find_dt();
      step(dt);
      t += dt;
    }
  }

  void step(double dt) {
    // extrapolate a half step in time using primitive equations
    primitive_update(0.5 * dt);

    // compute fluxes
    compute_fluxes();

    // update solultion from fluxes
    update_conserved(dt);

    // update primitive variables
    conserved_to_primitive();
  }

  void periodic_boundary(std::vector<double> &f) {
    // TODO: apply periodic boundary conditions to an array f
  }

  void primitive_update(double dt) {
    // TODO: update the primitive variables using Euler equations in primitive
    // form using an FTCS scheme
  }

  void extrapolate_to_interface() {
    // TODO: compute rho_L, rho_R, vx_L, vx_R, vy_L, vy_R, P_L, and P_R here
  }

  void compute_fluxes() {
    // TODO: compute the fluxes
  }

  void update_conserved(double dt) {
    // TODO: update the conserved variables using the fluxes
  }

  void output(int n) {
    std::ofstream outfile("output_rho_" + std::to_string(n) + ".csv");
    for (int j = 1; j < Ny - 1; j++) {
      for (int i = 1; i < Nx - 1; i++) {
        int idx = i + j * Nx;
        outfile << rho[idx];
        if (i != Nx - 2)
          outfile << ", ";
        else
          outfile << std::endl;
      }
    }
    outfile.close();
  }

  int Nx, Ny;
  double Lx, Ly;
  double dx, dy;
  double gamma, output_dt;
  std::vector<double> rho, vx, vy, P;                 // primitive variables
  std::vector<double> mass, mom_x, mom_y, energy;     // conserved variables
  // arrays to hold the results during primitive_update
  std::vector<double> rho_tmp, vx_tmp, vy_tmp, P_tmp;
  // arrays of fluxes for each conserved variable:
  std::vector<double> mass_flux_x, mass_flux_y;
  std::vector<double> momx_flux_x, momx_flux_y;
  std::vector<double> momy_flux_x, momy_flux_y;
  std::vector<double> energy_flux_x, energy_flux_y;
  // arrays for extrapolating to cell interfaces:
  std::vector<double> rho_Lx, rho_Ly, rho_Rx, rho_Ry;
  std::vector<double> vx_Lx, vx_Ly, vx_Rx, vx_Ry;
  std::vector<double> vy_Lx, vy_Ly, vy_Rx, vy_Ry;
  std::vector<double> P_Lx, P_Ly, P_Rx, P_Ry;
};
