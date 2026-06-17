#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <cmath>
#include <cassert>

using namespace std;

// ================= PARAMETERS =================
const int    N = 2000;

const double v0 = 2.0;
const double mu = 1.0;

const double sigma = 1.0;
const double r_int = sigma;

const double beta_CF  = 10;
const double alpha_CF = 20;
const double Dr = 0.0;

const double dt = 0.001;
const int    STEPS = 20000;
const int    SAVE_FREQ = 10;

// ================= GLOBALS ====================
double L;
double Rmax;
int    nCells;
double cell_size;

// ================= PARTICLE ==================
struct Particle {
    double x, y;
    double theta_n;
    double vx, vy;
};

// ================= PBC =======================
inline void wrap_position(double &x, double &y) {
    if (x < 0)  x += L;
    if (x >= L) x -= L;
    if (y < 0)  y += L;
    if (y >= L) y -= L;
}

inline void wrap_distance(double &dx, double &dy) {
    if (dx >  L/2) dx -= L;
    if (dx < -L/2) dx += L;
    if (dy >  L/2) dy -= L;
    if (dy < -L/2) dy += L;
}

// ================= SAVE ======================
void save_snapshot(ofstream &fout, const vector<Particle>& P, int frame) {
    fout << "FRAME " << frame << "\n";
    for (const auto &p : P) {
        fout << p.x << " " << p.y << " "
             << p.vx << " " << p.vy << " "
             << p.theta_n << "\n";
    }
}

// ================= INITIALIZE =================
vector<Particle> initialize_concentric_rings(int Ntarget, mt19937 &rng) {
    vector<Particle> P;
    P.reserve(Ntarget);

    uniform_real_distribution<double> ang(0.0, 2.0 * M_PI);

    int k = 0;
    while ((int)P.size() < Ntarget) {
        double Rk = (k + 0.5) * sigma;
        int Nk = max(1, int(2.0 * M_PI * Rk / sigma));

        for (int i = 0; i < Nk && (int)P.size() < Ntarget; i++) {
            double phi = 2.0 * M_PI * i / Nk;

            Particle p;
            p.x = Rk * cos(phi);
            p.y = Rk * sin(phi);

            // ---- orientation choice ----
            // p.theta_n = phi + M_PI/2.0;  // tangential
            // p.theta_n = phi;             // radial
            p.theta_n = ang(rng);          // random

            p.vx = v0 * cos(p.theta_n);
            p.vy = v0 * sin(p.theta_n);

            P.push_back(p);
        }
        k++;
    }

    Rmax = (k + 0.5) * sigma;
    return P;
}

// ================= CELL LIST =================
void build_cell_list(const vector<Particle>& P,
                     vector<vector<int>>& cell_list) {
    cell_list.assign(nCells * nCells, {});

    for (int i = 0; i < (int)P.size(); i++) {
        int cx = int(P[i].x / cell_size);
        int cy = int(P[i].y / cell_size);

        if (cx < 0) cx = 0;
        if (cy < 0) cy = 0;
        if (cx >= nCells) cx = nCells - 1;
        if (cy >= nCells) cy = nCells - 1;

        cell_list[cy*nCells + cx].push_back(i);
    }
}

// ================= UPDATE ====================
void update_step(vector<Particle>& P,
                 mt19937 &rng,
                 vector<vector<int>>& cell_list) {

    int Np = P.size();
    build_cell_list(P, cell_list);

    vector<pair<double,double>> forces(Np,{0,0});
    vector<pair<double,double>> JCF(Np,{0,0});

    normal_distribution<double> noise_dist(0.0,1.0);

    for (int cx = 0; cx < nCells; cx++) {
        for (int cy = 0; cy < nCells; cy++) {

            int cidx = cy*nCells + cx;
            const auto &plist = cell_list[cidx];

            for (int dxcell = -1; dxcell <= 1; dxcell++) {
                for (int dycell = -1; dycell <= 1; dycell++) {

                    int ncx = cx + dxcell;
                    int ncy = cy + dycell;
                    if (ncx < 0 || ncx >= nCells ||
                        ncy < 0 || ncy >= nCells)
                        continue;

                    int nidx = ncy*nCells + ncx;
                    if (nidx < cidx) continue;

                    const auto &nlist = cell_list[nidx];

                    for (int i : plist) {
                        for (int j : nlist) {
                            if (i >= j) continue;

                            double dx = P[j].x - P[i].x;
                            double dy = P[j].y - P[i].y;
                            wrap_distance(dx,dy);

                            double d2 = dx*dx + dy*dy;
                            if (d2 < 1e-12) continue;

                            double d = sqrt(d2);
                            if (d >= r_int) continue;

                            double dxn = dx / d;
                            double dyn = dy / d;

                            double Fmag = -beta_CF * (r_int/d - 1.0);
                            if (!isfinite(Fmag)) continue;

                            double Fx = Fmag * dxn;
                            double Fy = Fmag * dyn;

                            forces[i].first  += Fx;
                            forces[i].second += Fy;
                            forces[j].first  -= Fx;
                            forces[j].second -= Fy;

                            double qjx = cos(P[j].theta_n);
                            double qjy = sin(P[j].theta_n);
                            double w = 0.5 * (1.0 + dxn*qjx + dyn*qjy);

                            JCF[i].first  += w * dxn;
                            JCF[i].second += w * dyn;
                        }
                    }
                }
            }
        }
    }

    for (int i = 0; i < Np; i++) {
        P[i].vx = v0*cos(P[i].theta_n) + mu*forces[i].first;
        P[i].vy = v0*sin(P[i].theta_n) + mu*forces[i].second;
    }

    for (int i = 0; i < Np; i++) {
        double nx = cos(P[i].theta_n);
        double ny = sin(P[i].theta_n);
        double npx = -ny;
        double npy = nx;

        double noise = sqrt(2.0*Dr*dt) * noise_dist(rng);

        P[i].theta_n +=
            dt * alpha_CF * (JCF[i].first*npx + JCF[i].second*npy)
            + noise;
    }

    for (int i = 0; i < Np; i++) {
        P[i].x += P[i].vx * dt;
        P[i].y += P[i].vy * dt;
        wrap_position(P[i].x, P[i].y);
    }
}

// ================= MAIN ======================
int main() {
    mt19937 rng(12345);

    // temporary box
    L = 1.0;

    vector<Particle> P = initialize_concentric_rings(N, rng);

    // final box
    L = 4.0 * Rmax;

    // center configuration
    for (auto &p : P) {
        p.x += 0.5 * L;
        p.y += 0.5 * L;
    }

    cell_size = 2.5 * r_int;
    nCells = int(L / cell_size);
    if (nCells < 3)   nCells = 3;
    if (nCells > 200) nCells = 200;

    vector<vector<int>> cell_list;
    ofstream fout("simulation_data.txt");

    int frame = 0;
    for (int step = 0; step < STEPS; step++) {
        update_step(P, rng, cell_list);

        if (step % SAVE_FREQ == 0)
            save_snapshot(fout, P, frame++);
    }

    fout.close();
    cout << "Concentric-ring CF simulation finished.\n";
    return 0;
}
