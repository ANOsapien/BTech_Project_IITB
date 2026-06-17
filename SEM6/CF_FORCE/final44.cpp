#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <cmath>
using namespace std;
ofstream forceout;

// ================= PARAMETERS =================
const int    N        = 10000;
const double v0       = 1;
const double r_int    = 1.0;
const double beta_CF  = 1;
const double alpha_CF = 1;
const double Dr       = 0.0;
const double rho      = 1;
const double dt       = 0.001;
const int    STEPS    = 15000;
const int    SAVE_FREQ = 10;

// ================= GLOBALS ====================
double L;
int    nCells;
double cell_size;

// ================= PARTICLE ==================
struct Particle {
    double x, y;
    double theta_n;
    double vx, vy;
    double fx, fy;
};

ostream& operator<<(ostream& os, const Particle& p) {
    os << "Particle("
       << "x="       << p.x       << ", "
       << "y="       << p.y       << ", "
       << "theta_n=" << p.theta_n << ", "
       << "vx="      << p.vx      << ", "
       << "vy="      << p.vy
       << ")";
    return os;
}

// ================= PBC =======================
inline void wrap_position(double &x, double &y) {
    if (x <  0) x += L;
    if (x >= L) x -= L;
    if (y <  0) y += L;
    if (y >= L) y -= L;
}

// ================= I/O =======================
void save_snapshot(ofstream &fout, const vector<Particle>& P, int frame) {
    fout << "FRAME " << frame << "\n";
    for (const auto &p : P) {
        fout << p.x  << " " << p.y  << " "
             << p.vx << " " << p.vy << " "
             << p.theta_n << " "
             << p.fx << " " << p.fy << "\n";
    }
}

// ================= INITIALIZE =================
vector<Particle> initialize(int N, mt19937 &rng) {

    vector<Particle> P(N);

    int nSide = ceil(sqrt(N));      // lattice points per side
    double a = L / nSide;           // lattice spacing

    uniform_real_distribution<double> ang(0.0, 2*M_PI);

    int idx = 0;
    for (int i = 0; i < nSide && idx < N; i++) {
        for (int j = 0; j < nSide && idx < N; j++) {

            P[idx].x = (i + 0.5) * a;
            P[idx].y = (j + 0.5) * a;

            wrap_position(P[idx].x, P[idx].y);

            P[idx].theta_n = ang(rng);
            P[idx].vx = v0 * cos(P[idx].theta_n);
            P[idx].vy = v0 * sin(P[idx].theta_n);

            idx++;
        }
    }

    return P;
}


// ================= CELL LIST HELPERS ==========
inline int cell_id(int cx, int cy) {
    return cx + cy * nCells;
}

// ================= UPDATE (O(N) cell list) ====
void update(vector<Particle>& P, mt19937 &rng, int step) {
    int N = P.size();

    vector<pair<double,double>> forces(N, {0.0, 0.0});
    vector<pair<double,double>> JCF   (N, {0.0, 0.0});

    // ---------- build cell list ----------
    // Head-of-chain array + linked list: O(N), no inner allocations
    vector<int> head(nCells * nCells, -1);
    vector<int> next(N, -1);

    for (int i = 0; i < N; i++) {
        int cx = int(P[i].x / cell_size);
        int cy = int(P[i].y / cell_size);
        // clamp for floating-point edge cases
        if (cx >= nCells) cx = nCells - 1;
        if (cy >= nCells) cy = nCells - 1;
        int c   = cell_id(cx, cy);
        next[i] = head[c];
        head[c] = i;
    }

    // ---------- neighbor loop ----------
    // Only iterate over the 9-cell block once per ordered pair (i < j).
    // Strategy: for each cell c, pair all particles inside c with each other,
    // then pair them with particles in the 4 "forward" half-shell neighbors
    // (avoids double-counting while keeping PBC correct).
    //
    //   half-shell offsets relative to (cx, cy):
    //       (1,0), (1,1), (0,1), (-1,1)
    //   plus the cell itself for intra-cell pairs.

    const int dxc[] = { 0,  1,  1,  0, -1 };
    const int dyc[] = { 0,  0,  1,  1,  1 };
    // offset (0,0) handles intra-cell; the rest are the 4 forward neighbors.

    for (int cx = 0; cx < nCells; cx++) {
    for (int cy = 0; cy < nCells; cy++) {

        int c = cell_id(cx, cy);

        for (int k = 0; k < 5; k++) {

            int nx = (cx + dxc[k] + nCells) % nCells;
            int ny = (cy + dyc[k] + nCells) % nCells;
            int nc = cell_id(nx, ny);

            // iterate all pairs (i from cell c, j from cell nc)
            for (int i = head[c]; i != -1; i = next[i]) {
                // for intra-cell (k==0) start j from next[i] to avoid i==j
                // and double-counting; for cross-cell start from head[nc].
                int j_start = (k == 0) ? next[i] : head[nc];

                for (int j = j_start; j != -1; j = next[j]) {

                    // ------ PBC-correct signed displacement ------
                    double dx = P[j].x - P[i].x;
                    double dy = P[j].y - P[i].y;
                    if (dx >  L/2) dx -= L;
                    if (dx < -L/2) dx += L;
                    if (dy >  L/2) dy -= L;
                    if (dy < -L/2) dy += L;

                    double d2 = dx*dx + dy*dy;
                    if (d2 < 1e-12 || d2 >= r_int*r_int) continue;

                    double d   = sqrt(d2);
                    double dxn = dx / d;   // unit vector i → j
                    double dyn = dy / d;

                    // ------ contact repulsion (pushes i away from j) ------
                    double Fmag = beta_CF * (r_int / d - 1.0);
                    forces[i].first  -= Fmag * dxn;
                    forces[i].second -= Fmag * dyn;
                    forces[j].first  += Fmag * dxn;
                    forces[j].second += Fmag * dyn;

                    // ------ contact following (JCF) ------
                    // weight based on j's orientation relative to i→j direction
                    double qjx   = cos(P[j].theta_n);
                    double qjy   = sin(P[j].theta_n);
                    double weight = 0.5 * (1.0 + dxn*qjx + dyn*qjy);

                    JCF[i].first  += weight * dxn;
                    JCF[i].second += weight * dyn;
                }
            }
        }
    }}

    // ---------- velocity update + force log ----------
    for (int i = 0; i < N; i++) {
        const Particle& p = P[i];
        forceout << step         << ","
                 << p.x          << ","
                 << p.y          << ","
                 << p.vx         << ","
                 << p.vy         << ","
                 << p.theta_n    << ","
                 << forces[0].first  << ","
                 << forces[0].second << "\n";

        P[i].vx = v0 * cos(P[i].theta_n) + forces[i].first;
        P[i].vy = v0 * sin(P[i].theta_n) + forces[i].second;
        P[i].fx = forces[i].first;
        P[i].fy = forces[i].second;
    }

    // ---------- position update ----------
    for (int i = 0; i < N; i++) {
        P[i].x += P[i].vx * dt;
        P[i].y += P[i].vy * dt;
        wrap_position(P[i].x, P[i].y);
    }

    // ---------- angle update ----------
    for (int i = 0; i < N; i++) {
        double nx  =  cos(P[i].theta_n);
        double ny  =  sin(P[i].theta_n);
        double npx = -ny;   // perpendicular (left-hand normal)
        double npy =  nx;

        P[i].theta_n += dt * alpha_CF *
            (JCF[i].first * npx + JCF[i].second * npy);
    }
}

// ================= MAIN ======================
int main() {
    mt19937 rng(random_device{}());

    L = sqrt(N / rho);

    // cells must be at least r_int wide so 3x3 neighborhood covers all pairs
    cell_size = r_int;
    nCells    = max(1, int(L / cell_size));

    vector<Particle> P = initialize(N, rng);

    ofstream fout("simulation_data.txt");
    fout << "# L=" << L << "\n";

    forceout.open("particle1_full_trace.txt");
    forceout << "# step,x,y,vx,vy,theta,Fx,Fy\n";

    int frame = 0;
    for (int step = 0; step < STEPS; step++) {
        update(P, rng, step);
        if (step % SAVE_FREQ == 0)
            save_snapshot(fout, P, frame++);
    }

    fout.close();
    forceout.close();
    cout << "Simulation complete.\n";
    return 0;
}