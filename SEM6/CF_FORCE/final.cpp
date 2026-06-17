#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <cmath>
using namespace std;

// ================= PARAMETERS =================
const int    N = 20000;

const double v0 = 10.0;

const double r_int = 1.0;

const double beta_CF  = 400.0;
const double alpha_CF = 0.0;

const double Dr = 0.0;

const double rho = 0.5;

const double dt = 0.001;
const int    STEPS = 15000;
const int    SAVE_FREQ =10;

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
    for (const auto &p : P)
        fout << p.x << " " << p.y << " "
             << p.vx << " " << p.vy << " "
             << p.theta_n << "\n";
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

// ================= CELL LIST =================
void build_cell_list(const vector<Particle>& P,
                     vector<vector<int>>& cell_list) {

    cell_list.assign(nCells*nCells, {});

    for (int i=0;i<(int)P.size();i++) {
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

    for (int cx=0; cx<nCells; cx++) {
        for (int cy=0; cy<nCells; cy++) {

            const auto &plist = cell_list[cy*nCells+cx];

            for (int dxcell=-1; dxcell<=1; dxcell++) {
                for (int dycell=-1; dycell<=1; dycell++) {

                    int ncx=(cx+dxcell+nCells)%nCells;
                    int ncy=(cy+dycell+nCells)%nCells;
                    const auto &nlist = cell_list[ncy*nCells+ncx];

                    for (int i: plist) {
                        for (int j: nlist) {
                            if (i>=j) continue;

                            double dx=P[j].x-P[i].x;
                            double dy=P[j].y-P[i].y;
                            wrap_distance(dx,dy);

                            double d2=dx*dx+dy*dy;
                            if (d2<1e-12) continue;

                            double d=sqrt(d2);
                            if (d>=r_int) continue;

                            double dxn=dx/d;
                            double dyn=dy/d;

                            // --- nonlinear contact repulsion ---
                            double Fmag = -beta_CF*(r_int/d - 1.0);

                            double Fx = Fmag*dxn;
                            double Fy = Fmag*dyn;

                            forces[i].first  += Fx;
                            forces[i].second += Fy;
                            forces[j].first  -= Fx;
                            forces[j].second -= Fy;

                            // --- contact forcing (CF) ---
                            double qjx = cos(P[j].theta_n);
                            double qjy = sin(P[j].theta_n);

                            double weight = 0.5*(1.0 + dxn*qjx + dyn*qjy);

                            JCF[i].first  += weight*dxn;
                            JCF[i].second += weight*dyn;
                        }
                    }
                }
            }
        }
    }

    // velocity update
    for (int i=0;i<Np;i++) {
        P[i].vx = v0*cos(P[i].theta_n) + forces[i].first;
        P[i].vy = v0*sin(P[i].theta_n) + forces[i].second;
    }

    // angle update
    for (int i=0;i<Np;i++) {
        double nx=cos(P[i].theta_n);
        double ny=sin(P[i].theta_n);
        double npx=-ny;
        double npy= nx;

        P[i].theta_n += dt * alpha_CF *
            (JCF[i].first*npx + JCF[i].second*npy);
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

    mt19937 rng(random_device{}());

    // box size from density
    L = sqrt(N / rho);

    cell_size = 2.5 * r_int;
    nCells = max(3, int(L / cell_size));

    vector<Particle> P = initialize(N, rng);
    vector<vector<int>> cell_list;

    ofstream fout("simulation_data.txt");

    int frame=0;
    for (int step=0; step<STEPS; step++) {

        update_step(P, rng, cell_list);

        if (step % SAVE_FREQ == 0)
            save_snapshot(fout, P, frame++);
    }

    fout.close();

    cout << "Simulation complete.\n";
    return 0;
}
