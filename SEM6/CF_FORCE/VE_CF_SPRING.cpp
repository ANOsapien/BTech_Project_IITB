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

const double v0 = 1;
const double mu = 1.0;

const double sigma = 1.0;        // diameter
const double r_int = sigma;      // contact cutoff (center–center)

const double beta_CF = 1;      // nonlinear contact force strength
const double alpha_CF = 2;     // contact following strength
const double Dr = 0.0;          // angular noise

const double rho = 1;          // area fraction

const double dt = 0.01;
const int    STEPS = 20000;
const int    SAVE_FREQ = 10;

// ================= GLOBALS ====================
double L;
int nCells;
double cell_size;

// ================= PARTICLE ==================
struct Particle {
    double x, y;
    double theta_n;
    double vx, vy;
};

// ================= PBC =======================
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
// vector<Particle> initialize(int N, mt19937 &rng) {
//     vector<Particle> P(N);

//     int nSide = ceil(sqrt(N));
//     double dx = L / nSide;

//     uniform_real_distribution<double> jitter(-0.2*dx, 0.2*dx);
//     uniform_real_distribution<double> ang(0.0, 2*M_PI);

//     int idx = 0;
//     for (int i = 0; i < nSide && idx < N; i++) {
//         for (int j = 0; j < nSide && idx < N; j++) {

//             double x = (i + 0.5) * dx + jitter(rng);
//             double y = (j + 0.5) * dx + jitter(rng);

//             wrap_position(x, y);

//             P[idx].x = x;
//             P[idx].y = y;

//             P[idx].theta_n = ang(rng);
//             P[idx].vx = v0 * cos(P[idx].theta_n);
//             P[idx].vy = v0 * sin(P[idx].theta_n);

//             idx++;
//         }
//     }

//     return P;
// }


vector<Particle> initialize(int N, mt19937 &rng) {
    vector<Particle> P(N);

    double cx = 0.5 * L;
    double cy = 0.5 * L;

    // radius chosen so arc-length ~ sigma
    double R = (N * sigma) / (2.0 * M_PI);

    // optional small positional noise

    for (int i = 0; i < N; i++) {
        double phi = 2.0 * M_PI * i / N;

        P[i].x = cx + R * cos(phi) ;
        P[i].y = cy + R * sin(phi) ;

        wrap_position(P[i].x, P[i].y);

        // ---- orientation choice ----

        // (A) Tangential (circulating)
        P[i].theta_n = phi + M_PI / 2.0;

        // (B) Radial outward
        // P[i].theta_n = phi;

        // (C) Random
        // uniform_real_distribution<double> ang(0.0, 2*M_PI);
        // P[i].theta_n = ang(rng);

        P[i].vx = v0 * cos(P[i].theta_n);
        P[i].vy = v0 * sin(P[i].theta_n);
    }

    return P;
}


// ================= CELL LIST =================
void build_cell_list(const vector<Particle>& P,
                     vector<vector<int>>& cell_list) {
    int total = nCells * nCells;
    cell_list.assign(total, {});

    for (int i = 0; i < (int)P.size(); i++) {
        int cx = int(P[i].x / cell_size);
        int cy = int(P[i].y / cell_size);

        cx = (cx + nCells) % nCells;
        cy = (cy + nCells) % nCells;

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

                    int ncx = (cx + dxcell + nCells) % nCells;
                    int ncy = (cy + dycell + nCells) % nCells;
                    int nidx = ncy*nCells + ncx;
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

                            double invd = 1.0/d;
                            double dxn = dx*invd;
                            double dyn = dy*invd;

                            // ---- nonlinear contact force ----
                            double Fmag = -beta_CF * (r_int/d - 1.0);
                            double Fx = Fmag * dxn;
                            double Fy = Fmag * dyn;

                            forces[i].first  += Fx;
                            forces[i].second += Fy;
                            forces[j].first  -= Fx;
                            forces[j].second -= Fy;

                            // ---- contact following (CF) ----
                            double qjx = cos(P[j].theta_n);
                            double qjy = sin(P[j].theta_n);
                            double weight = 0.5 * (1.0 + dxn*qjx + dyn*qjy);

                            JCF[i].first  += weight * dxn;
                            JCF[i].second += weight * dyn;
                        }
                    }
                }
            }
        }
    }

    // velocity update
    for (int i=0;i<Np;i++) {
        P[i].vx = v0*cos(P[i].theta_n) + mu*forces[i].first;
        P[i].vy = v0*sin(P[i].theta_n) + mu*forces[i].second;
    }

    // angular update
    for (int i=0;i<Np;i++) {
        double nx = cos(P[i].theta_n);
        double ny = sin(P[i].theta_n);
        double npx = -ny;
        double npy = nx;

        double noise = sqrt(2.0*Dr*dt) * noise_dist(rng);

        P[i].theta_n +=
            dt * alpha_CF * (JCF[i].first*npx + JCF[i].second*npy)
            + noise;
    }

    // position update
    for (int i=0;i<Np;i++) {
        P[i].x += P[i].vx * dt;
        P[i].y += P[i].vy * dt;
        wrap_position(P[i].x,P[i].y);
    }
}

// ================= MAIN ======================
int main() {
    random_device rd;
    mt19937 rng(rd());

    // box size from density
    L = 4*R;

    cell_size = 2.5 * r_int;
    nCells = max(3, int(L / cell_size));

    vector<Particle> P = initialize(N, rng);
    vector<vector<int>> cell_list;

    ofstream fout("simulation_data.txt");

    int frame = 0;
    for (int step=0; step<STEPS; step++) {
        update_step(P, rng, cell_list);

        if (step % SAVE_FREQ == 0)
            save_snapshot(fout, P, frame++);
    }

    fout.close();
    cout << "CF + nonlinear contact simulation finished.\n";
    return 0;
}
