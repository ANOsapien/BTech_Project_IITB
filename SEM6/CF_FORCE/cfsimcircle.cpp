#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <cmath>
#include <cassert>

using namespace std;

// ================= PARAMETERS =================
const int    N = 200;

const double v0 = 2;
const double mu = 1.0;

const double sigma = 1.0;        // particle diameter
const double r_int = sigma;      // contact cutoff

const double beta_CF  = 1;     // nonlinear contact force
const double alpha_CF = 2;     // contact following strength
const double Dr = 0.0;           // angular noise

const double dt = 0.001;
const int    STEPS = 20000;
const int    SAVE_FREQ = 10;

// ================= GLOBALS ====================
double L;
double R;
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
vector<Particle> initialize_ring(int N) {
    vector<Particle> P(N);

    // double cx = 0.5 * L;
    // double cy = 0.5 * L;

    // sanity checks
    // assert(abs(2.0 * M_PI * R / N - sigma) < 1e-6);
    // assert(R < 0.5 * L);

    for (int i = 0; i < N; i++) {
        double phi = 2.0 * M_PI * i / N;

        P[i].x =  R * cos(phi);
        P[i].y =  R * sin(phi);

        // // tangential orientation
        P[i].theta_n = phi + M_PI / 2.0;
        // uniform_real_distribution<double> ang(0.0, 2*M_PI);
        // P[i].theta_n = ang(rng);
        // P[i].theta_n = phi;



        P[i].vx = v0 * cos(P[i].theta_n);
        P[i].vy = v0 * sin(P[i].theta_n);
    }

    return P;
}


// vector<Particle> initialize_ring(int N, mt19937 &rng) {
//     vector<Particle> P(N);

//     double cx = 0.5 * L;
//     double cy = 0.5 * L;

//     // exact touching-ring geometry checks
//     assert(abs(2.0 * M_PI * R / N - sigma) < 1e-6);
//     assert(R < 0.5 * L);

//     uniform_real_distribution<double> ang(0.0, 2.0 * M_PI);

//     for (int i = 0; i < N; i++) {
//         double phi = 2.0 * M_PI * i / N;

//         // position on ring
//         P[i].x = cx + R * cos(phi);
//         P[i].y = cy + R * sin(phi);

//         // random orientation
//         P[i].theta_n = ang(rng);

//         // velocity from orientation
//         P[i].vx = v0 * cos(P[i].theta_n);
//         P[i].vy = v0 * sin(P[i].theta_n);
//     }

//     return P;
// }


// ================= CELL LIST =================
void build_cell_list(const vector<Particle>& P,
                     vector<vector<int>>& cell_list) {
    int total = nCells * nCells;
    cell_list.assign(total, {});

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
                    if (nidx < cidx) continue;   // prevent double counting

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
                            if (d <= r_int) continue;

                            double dxn = dx / d;
                            double dyn = dy / d;

                            // ---- nonlinear contact force ----
                            double Fmag = -beta_CF * (r_int/d - 1.0);
                            if (!isfinite(Fmag)) continue;

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
    for (int i = 0; i < Np; i++) {
        P[i].vx = v0*cos(P[i].theta_n) + mu*forces[i].first;
        P[i].vy = v0*sin(P[i].theta_n) + mu*forces[i].second;
    }

    // angular update
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

    // position update
    for (int i = 0; i < Np; i++) {
        P[i].x += P[i].vx * dt;
        P[i].y += P[i].vy * dt;
        wrap_position(P[i].x, P[i].y);
    }
}

// ================= MAIN ======================
int main() {
    mt19937 rng(12345);

    // exact touching ring
    R = (N * sigma) / 2*(2.0 * M_PI);
    L = 10.0 * R;

    cell_size = 2.5 * r_int;
    nCells = int(L / cell_size);
    if (nCells < 3)   nCells = 3;
    if (nCells > 200) nCells = 200;

    vector<Particle> P = initialize_ring(N);
    vector<vector<int>> cell_list;

    ofstream fout("simulation_data.txt");

    int frame = 0;
    for (int step = 0; step < STEPS; step++) {
        update_step(P, rng, cell_list);

        if (step % SAVE_FREQ == 0)
            save_snapshot(fout, P, frame++);
    }

    fout.close();
    cout << "Single touching-ring CF simulation finished.\n";
    return 0;
}
