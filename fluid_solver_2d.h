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
    double vol = dx * dy;
    for (int idx = 0; idx < Nx * Ny; idx++) {
      int j = idx / Nx;
      double y = (j - 1) * dy;
      mass[idx] = rho[idx] * vol;
      mom_x[idx] = rho[idx] * vx[idx] * vol;
      mom_y[idx] = rho[idx] * vy[idx] * vol;
      energy[idx] = (P[idx] / (gamma - 1.0) +
                     rho[idx] * 0.5 * (vx[idx] * vx[idx] + vy[idx] * vy[idx])) *
                    vol;
    }
  }

  void conserved_to_primitive() {
    // TODO: Compute primitive variables from conserved ones
    double vol = dx * dy;
    for (int idx = 0; idx < Nx * Ny; idx++) {
      rho[idx] = mass[idx] / vol;
      vx[idx] = mom_x[idx] / mass[idx];
      vy[idx] = mom_y[idx] / mass[idx];
      P[idx] = (gamma - 1.0) *
               (energy[idx] / vol -
                0.5 * rho[idx] * (vx[idx] * vx[idx] + vy[idx] * vy[idx]));
    }
  }

  void init(const std::vector<double> &rho0, const std::vector<double> &vx0,
            const std::vector<double> &vy0, const std::vector<double> &P0) {
    // TODO: Initialize the primitive variables using the given initial
    // condition
    rho = rho0;
    vx = vx0;
    vy = vy0;
    P = P0;
    // obtain the conserved quantities from the primitive variables
    primitive_to_conserved();
  }

  double find_dt() {
    // TODO: Find the optimal dt that satisfies the CFL condition, and return
    // its value
    double dt = 100.0;
    for (int idx = 0; idx < Nx * Ny; idx++) {
      double dt_trial = 0.45 * std::min(dx, dy) /
                        (std::sqrt(gamma * P[idx] / rho[idx]) +
                         std::sqrt(vx[idx] * vx[idx] + vy[idx] * vy[idx]));
      if (dt_trial < dt)
        dt = dt_trial;
    }
    return dt;
  }

  void solve(double t0, double t_end) {
    // Solve the fluid equations, starting from t0 and stoping at t_end
    double t = t0;
    int n = 0; // n labels the output file
    while (t < t_end) {
      std::cout << "t = " << t << ", n = " << n << std::endl;
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

    // extrapolate primitive variables to interfaces
    extrapolate_to_interface();

    // compute fluxes
    compute_fluxes();

    // update solultion from fluxes
    update_conserved(dt);

    // update primitive variables
    conserved_to_primitive();
  }

  void periodic_boundary(std::vector<double> &f) {
    // TODO: apply periodic boundary conditions to an array f
    for (int i = 1; i < Nx - 1; i++) {
      f[i] = f[i + (Ny - 2) * Nx];
      f[i + (Ny - 1) * Nx] = f[i + Nx];
    }
    for (int j = 1; j < Ny - 1; j++) {
      f[j * Nx] = f[j * Nx + Nx - 2];
      f[j * Nx + Nx - 1] = f[j * Nx + 1];
    }
  }

  void primitive_update(double dt) {
    // TODO: update the primitive variables using Euler equations in primitive
    // form using an FTCS scheme
    for (int j = 1; j < Ny - 1; j++) {
      for (int i = 1; i < Nx - 1; i++) {
        int idx = i + j * Nx;
        rho_tmp[idx] =
            rho[idx] -
            dt * (vx[idx] * (rho[idx + 1] - rho[idx - 1]) / (2.0 * dx) +
                  vy[idx] * (rho[idx + Nx] - rho[idx - Nx]) / (2.0 * dy) +
                  rho[idx] * ((vx[idx + 1] - vx[idx - 1]) / (2.0 * dx) +
                              (vy[idx + Nx] - vy[idx - Nx]) / (2.0 * dy)));
        vx_tmp[idx] =
            vx[idx] -
            dt * (vx[idx] * (vx[idx + 1] - vx[idx - 1]) / (2.0 * dx) +
                  vy[idx] * (vx[idx + Nx] - vx[idx - Nx]) / (2.0 * dy) +
                  (P[idx + 1] - P[idx - 1]) / (2.0 * dx) / rho[idx]);

        vy_tmp[idx] =
            vy[idx] -
            dt * (vx[idx] * (vy[idx + 1] - vy[idx - 1]) / (2.0 * dx) +
                  vy[idx] * (vy[idx + Nx] - vy[idx - Nx]) / (2.0 * dy) +
                  (P[idx + Nx] - P[idx - Nx]) / (2.0 * dy) / rho[idx]);

        P_tmp[idx] =
            P[idx] - dt * (vx[idx] * (P[idx + 1] - P[idx - 1]) / (2.0 * dx) +
                           vy[idx] * (P[idx + Nx] - P[idx - Nx]) / (2.0 * dy) +
                           gamma * P[idx] *
                               ((vx[idx + 1] - vx[idx - 1]) / (2.0 * dx) +
                                (vy[idx + Nx] - vy[idx - Nx]) / (2.0 * dy)));
      }
    }

    periodic_boundary(rho_tmp);
    periodic_boundary(vx_tmp);
    periodic_boundary(vy_tmp);
    periodic_boundary(P_tmp);
    rho = rho_tmp;
    vx = vx_tmp;
    vy = vy_tmp;
    P = P_tmp;
  }

  void extrapolate_to_interface() {
    // TODO: compute rho_L, rho_R, vx_L, vx_R, vy_L, vy_R, P_L, and P_R here
    for (int j = 1; j < Ny - 1; j++) {
      for (int i = 1; i < Nx - 1; i++) {
        int idx = i + j * Nx;
        // L means left face, R means right face
        rho_Lx[idx] = rho_tmp[idx] - (rho[idx + 1] - rho[idx - 1]) / 4.0;
        rho_Rx[idx] = rho_tmp[idx] + (rho[idx + 1] - rho[idx - 1]) / 4.0;
        vx_Lx[idx] = vx_tmp[idx] - (vx[idx + 1] - vx[idx - 1]) / 4.0;
        vx_Rx[idx] = vx_tmp[idx] + (vx[idx + 1] - vx[idx - 1]) / 4.0;
        vy_Lx[idx] = vy_tmp[idx] - (vy[idx + 1] - vy[idx - 1]) / 4.0;
        vy_Rx[idx] = vy_tmp[idx] + (vy[idx + 1] - vy[idx - 1]) / 4.0;
        P_Lx[idx] = P_tmp[idx] - (P[idx + 1] - P[idx - 1]) / 4.0;
        P_Rx[idx] = P_tmp[idx] + (P[idx + 1] - P[idx - 1]) / 4.0;
        rho_Ly[idx] = rho_tmp[idx] - (rho[idx + Nx] - rho[idx - Nx]) / 4.0;
        rho_Ry[idx] = rho_tmp[idx] + (rho[idx + Nx] - rho[idx - Nx]) / 4.0;
        vx_Ly[idx] = vx_tmp[idx] - (vx[idx + Nx] - vx[idx - Nx]) / 4.0;
        vx_Ry[idx] = vx_tmp[idx] + (vx[idx + Nx] - vx[idx - Nx]) / 4.0;
        vy_Ly[idx] = vy_tmp[idx] - (vy[idx + Nx] - vy[idx - Nx]) / 4.0;
        vy_Ry[idx] = vy_tmp[idx] + (vy[idx + Nx] - vy[idx - Nx]) / 4.0;
        P_Ly[idx] = P_tmp[idx] - (P[idx + Nx] - P[idx - Nx]) / 4.0;
        P_Ry[idx] = P_tmp[idx] + (P[idx + Nx] - P[idx - Nx]) / 4.0;
      }
    }
    periodic_boundary(rho_Lx);
    periodic_boundary(rho_Rx);
    periodic_boundary(rho_Ly);
    periodic_boundary(rho_Ry);
    periodic_boundary(vx_Lx);
    periodic_boundary(vx_Rx);
    periodic_boundary(vx_Ly);
    periodic_boundary(vx_Ry);
    periodic_boundary(vy_Lx);
    periodic_boundary(vy_Rx);
    periodic_boundary(vy_Ly);
    periodic_boundary(vy_Ry);
    periodic_boundary(P_Lx);
    periodic_boundary(P_Rx);
    periodic_boundary(P_Ly);
    periodic_boundary(P_Ry);
  }

  void compute_fluxes() {
    // TODO: compute the fluxes
    for (int j = 1; j < Ny - 1; j++) {
      for (int i = 1; i < Nx - 1; i++) {
        int idx = i + j * Nx;

        // compute x fluxes first
        double c_l = std::sqrt(gamma * P_Lx[idx + 1] / rho_Lx[idx + 1]) +
                     std::abs(vx_Lx[idx + 1]);
        double c_r =
            std::sqrt(gamma * P_Rx[idx] / rho_Rx[idx]) + std::abs(vx_Rx[idx]);
        double c = std::max(c_l, c_r);

        double en_L = P_Lx[idx + 1] / (gamma - 1.0) +
                      0.5 * rho_Lx[idx + 1] *
                          (square(vx_Lx[idx + 1]) + square(vy_Lx[idx + 1]));
        double en_R =
            P_Rx[idx] / (gamma - 1.0) +
            0.5 * rho_Rx[idx] * (square(vx_Rx[idx]) + square(vy_Rx[idx]));

        // double P_mid = P_Lx[idx + 1]

        mass_flux_x[idx] =
            0.5 * (rho_Rx[idx] * vx_Rx[idx] + rho_Lx[idx + 1] * vx_Lx[idx + 1]) * dy +
            0.5 * c * (rho_Rx[idx] - rho_Lx[idx + 1]) * dy;
        momx_flux_x[idx] =
            // (momx_mid * momx_mid / rho_mid + P_mid) * dy +
            0.5 *
                (rho_Rx[idx] * vx_Rx[idx] * vx_Rx[idx] + P_Rx[idx] +
                 rho_Lx[idx + 1] * vx_Lx[idx + 1] * vx_Lx[idx + 1] + P_Lx[idx + 1]) *
                dy +
            0.5 * c *
                (rho_Rx[idx] * vx_Rx[idx] - rho_Lx[idx + 1] * vx_Lx[idx + 1]) * dy;
        momy_flux_x[idx] =
            // (momx_mid * momy_mid / rho_mid) * dy +
            0.5 *
                (rho_Rx[idx] * vx_Rx[idx] * vy_Rx[idx] +
                 rho_Lx[idx + 1] * vx_Lx[idx + 1] * vy_Lx[idx + 1]) *
                dy +
            0.5 * c *
                (rho_Rx[idx] * vy_Rx[idx] - rho_Lx[idx + 1] * vy_Lx[idx + 1]) * dy;
        // energy_flux_x[idx] = ((en_mid + P_mid) * momx_mid / rho_mid) * dy -
        energy_flux_x[idx] = 0.5 * ((en_L + P_Lx[idx + 1]) * vx_Lx[idx + 1] +
                                    (en_R + P_Rx[idx]) * vx_Rx[idx]) * dy -
                             0.5 * c * (en_L - en_R) * dy;

        // flux in y
        c_r = std::sqrt(gamma * P_Ly[idx + Nx] / rho_Ly[idx + Nx]) +
                     std::abs(vy_Ly[idx + Nx]);
        c_l =
            std::sqrt(gamma * P_Ry[idx] / rho_Ry[idx]) + std::abs(vy_Ry[idx]);
        c = std::max(c_l, c_r);

        en_L = P_Ly[idx + Nx] / (gamma - 1.0) +
                      0.5 * rho_Ly[idx + Nx] *
                          (vx_Ly[idx + Nx] * vx_Ly[idx + Nx] +
                           vy_Ly[idx + Nx] * vy_Ly[idx + Nx]);
        en_R =
            P_Ry[idx] / (gamma - 1.0) +
            0.5 * rho_Ry[idx] * (vx_Ry[idx] * vx_Ry[idx] + vy_Ry[idx] * vy_Ry[idx]);

        mass_flux_y[idx] =
            0.5 * (rho_Ry[idx] * vy_Ry[idx] + rho_Ly[idx + Nx] * vy_Ly[idx + Nx]) * dx +
            0.5 * c * (rho_Ry[idx] - rho_Ly[idx + Nx]) * dx;
        momx_flux_y[idx] =
            // (momy_mid * momx_mid / rho_mid) * dx +
            0.5 *
                (rho_Ry[idx] * vx_Ry[idx] * vy_Ry[idx] +
                 rho_Ly[idx + Nx] * vx_Ly[idx + Nx] * vy_Ly[idx + Nx]) *
                dx +
            0.5 * c *
                (rho_Ry[idx] * vx_Ry[idx] - rho_Ly[idx + Nx] * vx_Ly[idx + Nx]) *
                dx;
        momy_flux_y[idx] =
            // (momy_mid * momy_mid / rho_mid + P_mid) * dx +
            0.5 *
                (rho_Ry[idx] * vy_Ry[idx] * vy_Ry[idx] + P_Ry[idx] +
                 rho_Ly[idx + Nx] * vy_Ly[idx + Nx] * vy_Ly[idx + Nx] + P_Ly[idx + Nx]) *
                dx +
            0.5 * c *
                (rho_Ry[idx] * vy_Ry[idx] - rho_Ly[idx + Nx] * vy_Ly[idx + Nx]) *
                dx;
        energy_flux_y[idx] = 0.5 * ((en_L + P_Ly[idx + Nx]) * vy_Ly[idx + Nx] +
                                    (en_R + P_Ry[idx]) * vy_Ry[idx]) * dx -
                             0.5 * c * (en_L - en_R) * dx;

      }
    }

    periodic_boundary(mass_flux_x);
    periodic_boundary(mass_flux_y);
    periodic_boundary(momx_flux_x);
    periodic_boundary(momx_flux_y);
    periodic_boundary(momy_flux_x);
    periodic_boundary(momy_flux_y);
    periodic_boundary(energy_flux_x);
    periodic_boundary(energy_flux_y);
  }

  void update_conserved(double dt) {
    // TODO: update the conserved variables using the fluxes
    for (int j = 1; j < Ny - 1; j++) {
      for (int i = 1; i < Nx - 1; i++) {
        int idx = i + j * Nx;
        mass[idx] -= dt * (mass_flux_x[idx] - mass_flux_x[idx - 1]) +
                     dt * (mass_flux_y[idx] - mass_flux_y[idx - Nx]);
        mom_x[idx] -= dt * (momx_flux_x[idx] - momx_flux_x[idx - 1]) +
                      dt * (momx_flux_y[idx] - momx_flux_y[idx - Nx]);
        mom_y[idx] -= dt * (momy_flux_x[idx] - momy_flux_x[idx - 1]) +
                      dt * (momy_flux_y[idx] - momy_flux_y[idx - Nx]);
        energy[idx] -= dt * (energy_flux_x[idx] - energy_flux_x[idx - 1]) +
                       dt * (energy_flux_y[idx] - energy_flux_y[idx - Nx]);
      }
    }
    periodic_boundary(mass);
    periodic_boundary(mom_x);
    periodic_boundary(mom_y);
    periodic_boundary(energy);
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
