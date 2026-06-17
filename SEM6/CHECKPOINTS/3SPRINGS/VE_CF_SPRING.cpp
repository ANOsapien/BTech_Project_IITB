#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <cassert>
#include <string>
#include <sstream>
#include <filesystem>
using namespace std;


// =============== PARAMETERS ===============
const double v0   = 1;
const double alpha_CF=0;
const double tau  = 1e12;

const double mu   = 1.0;
const double sigma = 1;
const double Dr = 0.0;
const double k = 10.0;
const int N = 3;
const double rho = 0.1;
const double radius= sigma/2;


//Pe=v0/(sigma*Dr)
//cf=alpha_CF/Dr
//g=1/(tau*Dr)

// varying params are just: v0, alpha_CF & tau
const double pi=22.0/7.0;
const double dt = 0.001;
const int STEPS = 10000;
const int BURN_IN_STEPS = 1;
const int SAVE_FREQ = 10; // Save every N steps after burn-in

// =============== GLOBALS for neighbor list ===============
double L;
int nCells;
double cell_size;

// =============== PARTICLE STRUCT ===============
struct Particle {
    double x, y;
    double theta_n;
    double vx, vy;
};

// =============== UTILITIES ===============
inline void wrap_position(double &x, double &y) {
    if (x < 0) x += L;
    if (x >= L) x -= L;
    if (y < 0) y += L;
    if (y >= L) y -= L;
}


inline void wrap_distance(double &dx, double &dy) {
    if (dx >  L/2) dx -= L;
    if (dx < -L/2) dx += L;
    if (dy >  L/2) dy -= L;
    if (dy < -L/2) dy += L;
}

inline int mod_wrap(int idx, int m) {
    int r = idx % m;
    if (r < 0) r += m;
    return r;
}

// =============== SAVE SNAPSHOT ===============
void save_snapshot(ofstream &fout, const vector<Particle>& P, int frame) {
    fout << "FRAME " << frame << "\n";
    for (const auto &p : P) {
        fout << p.x << " " << p.y << " " << p.vx << " " 
             << p.vy << " " << p.theta_n << "\n";
    }
    fout.flush();
}

//=============== INITIALIZE ===============

// vector<Particle> initialize(int N, mt19937 &rng) {
//     int n = int(sqrt(N));
//     double alat = 2 * L / (2 * n + 3);
//     uniform_real_distribution<double> ang(0.0, 2 * M_PI);

    
//     vector<Particle> P(N);
//     int nptl = 0;
    
//     for (int i = 0; i < n; i++) {
//         for (int j = 0; j < n; j++) {
//             if (j % 2 == 0)
//                 P[nptl].x = L/2.0 + i * alat;  // Adjusted for [0,L)
//             else
//                 P[nptl].x = L/2.0 + (2.0*i + 1.0)*alat/2.0;
            
//             if (i == 0)
//                 P[nptl].y = L/2.0 + (2.0*j + 1.0)*alat/2.0;
//             else
//                 P[nptl].y = L/2.0 + j * alat;
            
//             P[nptl].theta_n =ang(rng);
//             P[nptl].vx = v0 * cos(P[nptl].theta_n);
//             P[nptl].vy = v0 * sin(P[nptl].theta_n);
//             nptl++;
//         }
//     }
//     return P;
//  }




// vector<Particle> initialize(int N, mt19937 &rng) {
//     int nSide = ceil(sqrt(N));
//     double dx = L / nSide;

//     uniform_real_distribution<double> ang(0.0, 2 * M_PI);
//     vector<Particle> P(N);

//     int idx = 0;
//     for (int i = 0; i < nSide && idx < N; i++) {
//         for (int j = 0; j < nSide && idx < N; j++) {
//             P[idx].x = (i + 0.5) * dx;
//             P[idx].y = (j + 0.5) * dx;
//             P[idx].theta_n = ang(rng);
//             P[idx].vx = v0 * cos(P[idx].theta_n);
//             P[idx].vy = v0 * sin(P[idx].theta_n);
//             idx++;
//         }
//     }
//     return P;
// }


vector<Particle> initialize(int N, mt19937 &rng) {
    vector<Particle> P(3);

    double a = sigma;     // spring rest length
    L = 10.0 * a;               // large box to avoid PBC

    double cx = L / 2.0;
    double cy = L / 2.0;

    // Equilateral triangle
    P[0].x = cx - a/2;
    P[0].y = cy - sqrt(3)*a/6;

    P[1].x = cx + a/2;
    P[1].y = cy - sqrt(3)*a/6;

    P[2].x = cx;
    P[2].y = cy + sqrt(3)*a/3;

    // Orient along bonds (cyclic)
    auto angle_to = [&](int i, int j){
        return atan2(P[j].y - P[i].y, P[j].x - P[i].x);
    };

    P[0].theta_n = angle_to(0,1);
    P[1].theta_n = angle_to(1,2);
    P[2].theta_n = angle_to(2,0);

    for (int i = 0; i < 3; i++) {
        P[i].vx = v0 * cos(P[i].theta_n);
        P[i].vy = v0 * sin(P[i].theta_n);
    }

    return P;
}


// =============== NEIGHBOR LIST ===============
void build_cell_list(const vector<Particle>& P, vector<vector<int>>& cell_list) {
    int total = nCells * nCells;
    cell_list.assign(total, {});
    
    int avg = max(1, int(P.size() / (total)));
    for (int i=0;i<total;i++) cell_list[i].reserve(avg+1);

    for (int i=0;i<(int)P.size(); ++i) {
        int cx = int(floor(P[i].x / cell_size));
        int cy = int(floor(P[i].y / cell_size));
        if (cx < 0) cx = 0;
        if (cx >= nCells) cx = nCells - 1;
        if (cy < 0) cy = 0;
        if (cy >= nCells) cy = nCells - 1;
        int idx = cy * nCells + cx;
        assert(idx >= 0 && idx < total);
        cell_list[idx].push_back(i);
    }
}

// =============== UPDATE STEP ===============

void update_step(vector<Particle>& P, mt19937 &rng, vector<vector<int>>& cell_list) {
    int N = (int)P.size();
    build_cell_list(P, cell_list);

    // === Force and alignment accumulators ===
    static vector<pair<double,double>> forces(N, {0.0, 0.0});
    static vector<pair<double,double>> JCF(N, {0.0, 0.0});
    static vector<double> align_sum(N, 0.0);

    for (int i = 0; i < N; ++i) {
        forces[i] = {0.0, 0.0};
        JCF[i] = {0.0, 0.0};
        align_sum[i] = 0.0;
    }

    // === Compute interactions using neighbor list ===
    for (int cx = 0; cx < nCells; ++cx) {
        for (int cy = 0; cy < nCells; ++cy) {
            int cindex = cy * nCells + cx;
            const auto &plist = cell_list[cindex];

            for (int ddx = -1; ddx <= 1; ++ddx) {
                int ncx = mod_wrap(cx + ddx, nCells);
                for (int ddy = -1; ddy <= 1; ++ddy) {
                    int ncy = mod_wrap(cy + ddy, nCells);
                    int nindex = ncy * nCells + ncx;
                    const auto &nlist = cell_list[nindex];

                    for (int i : plist) {
                        for (int j : nlist) {
                            if (i >= j) continue;

                            double dx = P[j].x - P[i].x;
                            double dy = P[j].y - P[i].y;
                            wrap_distance(dx, dy);
                            double d2 = dx*dx + dy*dy;
                            if (d2 < 1e-24) continue;
                            double d = sqrt(d2);
                            double invd = 1.0 / d;

                            // === Spring-like contact force (for r < 3σ) ===
                            // if (d < 2*sigma) {
                            double Fmag = -k * (2*sigma - d);
                            double Fx = Fmag * dx * invd;
                             double Fy = Fmag * dy * invd;
                            forces[i].first += Fx;
                            forces[i].second += Fy;
                            forces[j].first -= Fx;
                            forces[j].second -= Fy;
                            

                            // === Angular alignment (for r < 2σ) ===
                            if (d < sigma) {
                                align_sum[i] += sin(P[j].theta_n - P[i].theta_n);
                                align_sum[j] += sin(P[i].theta_n - P[j].theta_n);
                            }

                            // === Contact forcing (as before) ===
                            if (d <= 2*sigma) {
                                double dxn = dx * invd;
                                double dyn = dy * invd;
                                double qjx = cos(P[j].theta_n);
                                double qjy = sin(P[j].theta_n);
                                double dot = dxn * qjx + dyn * qjy;
                                double weight = 0.5 * (1.0 + dot);
                                JCF[i].first += weight * dxn;
                                JCF[i].second += weight * dyn;
                            }
                        }
                    }
                }
            }
        }
    }

    // === Velocity update ===
    for (int i = 0; i < N; ++i) {
        double vx_prop = (v0 ) * cos(P[i].theta_n);
        double vy_prop = (v0 ) * sin(P[i].theta_n);
        P[i].vx = vx_prop + (mu ) * forces[i].first;
        P[i].vy = vy_prop + (mu ) * forces[i].second;
    }

    // === Angular update: alignment + contact forcing + noise ===
    normal_distribution<double> noise_dist(0.0, 1.0);
    vector<double> new_theta(N);

    for (int i = 0; i < N; ++i) {
        // Contact force torque
        double nx = cos(P[i].theta_n);
        double ny = sin(P[i].theta_n);
        double nperp_x = -ny;
        double nperp_y = nx;
        double torque_CF = JCF[i].first * nperp_x + JCF[i].second * nperp_y;

        // Alignment term + contact torque + Gaussian noise
        double noise = sqrt(2.0 * Dr * dt) * noise_dist(rng);
        new_theta[i] = P[i].theta_n
            + (dt/tau)* align_sum[i]
            + (dt * (alpha_CF) * torque_CF)
            + noise;
    }

    // === Position and angle updates ===
    for (int i = 0; i < N; ++i) {
        P[i].theta_n = new_theta[i];
        P[i].x += P[i].vx * dt;
        P[i].y += P[i].vy * dt;
        wrap_position(P[i].x, P[i].y);
    }
}


// =============== RUN SIMULATION ===============
void run_simulation(int N, double rho, mt19937 &rng) {
    L = sqrt(N * pi * radius * radius / rho);    
    cell_size = 3*sigma;
    nCells = max(3, int(floor(L / cell_size)));

    vector<Particle> P = initialize(N, rng);
    vector<vector<int>> cell_list;

    ofstream fout("simulation_data.txt");
    fout << "# N=" << N << " rho=" << rho << " L=" << L << "\n";

    int frame = 0;
    for (int step = 0; step < STEPS; ++step) {
        update_step(P, rng, cell_list);
        if (step >= BURN_IN_STEPS && (step - BURN_IN_STEPS) % SAVE_FREQ == 0) {
            save_snapshot(fout, P, frame++);
        }
    }
    fout.close();
    cout << "Simulation complete. Saved " << frame << " frames to simulation_data.txt\n";
}

// =============== MAIN ===============
int main() {
    random_device rd;
    mt19937 rng(rd());
    run_simulation(N, rho, rng);
    return 0;
}