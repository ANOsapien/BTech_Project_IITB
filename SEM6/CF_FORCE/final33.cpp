#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <cmath>
using namespace std;
ofstream forceout;
// ================= PARAMETERS =================
const int    N =2;

const double v0 = 1;

const double r_int = 1.0;

const double beta_CF  = 1;
const double alpha_CF = 0;

const double Dr = 0.0;

const double rho = 0.01;

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
    double fx, fy;
};

ostream& operator<<(ostream& os, const Particle& p) {
    os << "Particle("
       << "x=" << p.x << ", "
       << "y=" << p.y << ", "
       << "theta_n=" << p.theta_n << ", "
       << "vx=" << p.vx << ", "
       << "vy=" << p.vy
       << ")";
    return os;
}

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


void save_snapshot(ofstream &fout, const vector<Particle>& P, int frame) {

    fout << "FRAME " << frame << "\n";

    // columns: x y vx vy theta
    for (const auto &p : P) {
        fout << p.x << " "
             << p.y << " "
             << p.vx << " "
             << p.vy << " "
             << p.theta_n << " "
             << p.fx << " "
             << p.fy <<  "\n";
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


vector<Particle> initialize(int N) {

    vector<Particle> P(N);

    double R = 0.35 * L;          // ring radius
    double cx = 0.5 * L;
    double cy = 0.5 * L;

    for (int i = 0; i < N; i++) {

        double phi = 2.0 * M_PI * i / N;

        P[i].x = cx + R * cos(phi);
        P[i].y = cy + R * sin(phi);
        wrap_position(P[i].x, P[i].y);

        // tangential direction
        P[i].theta_n = phi ;

        P[i].vx = -v0 * cos(P[i].theta_n);
        P[i].vy = -v0 * sin(P[i].theta_n);
    }

    return P;
}


// vector<Particle> initialize(int N) {

//     vector<Particle> P(N);

//     double R  = 0.35 * L;
//     double cx = 0.5 * L;
//     double cy = 0.5 * L;

//     double alpha = 0.4;   // spiral strength (0–1 typical)

//     for (int i = 0; i < N; i++) {

//         double phi = 2.0 * M_PI * i / N;

//         // position on ring
//         P[i].x = cx + R * cos(phi);
//         P[i].y = cy + R * sin(phi);
//         wrap_position(P[i].x, P[i].y);

//         // spiral direction = tangential + radial
//         double theta_tan = phi + M_PI/2.0;
//         double theta_rad = phi;

//         double vx = cos(theta_tan) + alpha * cos(theta_rad);
//         double vy = sin(theta_tan) + alpha * sin(theta_rad);

//         // normalize to speed v0
//         double norm = sqrt(vx*vx + vy*vy);
//         P[i].vx = v0 * vx / norm;
//         P[i].vy = v0 * vy / norm;

//         P[i].theta_n = atan2(P[i].vy, P[i].vx);
//     }

//     return P;
// }

// vector<Particle> initialize(int N, mt19937 &rng) {
//     vector<Particle> P(N);

//     uniform_real_distribution<double> pos(0.0, L);
//     uniform_real_distribution<double> ang(0.0, 2.0 * M_PI);

//     for (int i = 0; i < N; i++) {
//         P[i].x = pos(rng);
//         P[i].y = pos(rng);

//         P[i].theta_n = ang(rng);

//         P[i].vx = v0 * cos(P[i].theta_n);
//         P[i].vy = v0 * sin(P[i].theta_n);
//     }

//     return P;
// }

bool nearly_equal(double a, double b,
                  double abs_tol = 1e-12,
                  double rel_tol = 1e-9)
{
    return std::fabs(a - b) <= std::max(abs_tol, rel_tol * std::max(std::fabs(a), std::fabs(b)));
}


void update(vector<Particle>& P, mt19937 &rng, int step){
    int N = P.size();
    vector<pair<double,double>> forces(N);
    vector<pair<double,double>> JCF(N);
    // forces[i].first = x directions
    for(int i=0; i<N; i++){
        for(int j=i+1; j<N; j++){
            
double dx = P[j].x - P[i].x;   // signed
double dy = P[j].y - P[i].y;   // signed

// PBC: wrap to [-L/2, L/2]
if (dx >  L/2) dx -= L;
if (dx < -L/2) dx += L;
if (dy >  L/2) dy -= L;
if (dy < -L/2) dy += L;

double d2 = dx*dx + dy*dy;
if (d2 < 1e-12 || d2 >= r_int*r_int) continue;

double d   = sqrt(d2);
double dxn = dx / d;   // now correctly signed: points from i toward j
double dyn = dy / d;

// Repulsion: push i away from j  →  force on i is in -rij direction
double Fmag = beta_CF * (r_int/d - 1.0);
forces[i].first  -= Fmag * dxn;   // note the minus
forces[i].second -= Fmag * dyn;
forces[j].first  += Fmag * dxn;
forces[j].second += Fmag * dyn;

// JCF: dxn/dyn now point in the correct direction
double qjx   = cos(P[j].theta_n);
double qjy   = sin(P[j].theta_n);
double weight = 0.5*(1.0 + dxn*qjx + dyn*qjy);  // dxn points i→j
JCF[i].first  += weight * dxn;
JCF[i].second += weight * dyn;
            
        }
    }

    for (int i = 0; i < N; i++) {
        const Particle& p = P[i];

        forceout << step << ","
                << p.x << ","
                << p.y << ","
                << p.vx << ","
                << p.vy << ","
                << p.theta_n << ","
                << forces[0].first << ","
                << forces[0].second << "\n";
        P[i].vx = v0*cos(P[i].theta_n) + forces[i].first;
        P[i].vy = v0*sin(P[i].theta_n) + forces[i].second;
        P[i].fx = forces[i].first;
        P[i].fy = forces[i].second;
    }
    float tfx = 0, tfy = 0;
    // ---- position update ----
    for (int i = 0; i < N; i++) {
        P[i].x += P[i].vx * dt;
        P[i].y += P[i].vy * dt;
        wrap_position(P[i].x, P[i].y);
        tfx += P[i].fx;
        tfy += P[i].fy;
    }
    
   // cout << tfx << ' ' << tfy << endl;

    // ---- angle update ----
    for (int i = 0; i < N; i++) {
        double nx = cos(P[i].theta_n);
        double ny = sin(P[i].theta_n);
        double npx = -ny;
        double npy =  nx;

        P[i].theta_n += dt * alpha_CF *
            (JCF[i].first*npx + JCF[i].second*npy);
    }
    if(nearly_equal(tfx,0) || nearly_equal(tfy,0)){
        cout << "wronggggg at step number " << step << endl;
    }
    // if(step%100==0) cout << forces[0].first << ' ' << forces[0].second << '\n';
    
}


// inline int cell_id(int cx, int cy) {
//     return cx + cy * nCells;
// }
// void update(vector<Particle>& P, mt19937 &rng, int step){

//     int N = P.size();

//     vector<pair<double,double>> forces(N, {0.0,0.0});
//     vector<pair<double,double>> JCF(N, {0.0,0.0});

//     // ---------------- build cell list ----------------
//     vector<vector<int>> cell_list(nCells * nCells);

//     for (int i = 0; i < N; i++) {
//         int cx = int(P[i].x / cell_size);
//         int cy = int(P[i].y / cell_size);

//         if (cx >= nCells) cx = nCells - 1;
//         if (cy >= nCells) cy = nCells - 1;

//         cell_list[cell_id(cx, cy)].push_back(i);
//     }

//     // ---------------- neighbor loop (O(N)) ----------------
//     for (int cx = 0; cx < nCells; cx++) {
//     for (int cy = 0; cy < nCells; cy++) {

//         int c = cell_id(cx, cy);

//         for (int i : cell_list[c]) {

//             for (int dx_cell = -1; dx_cell <= 1; dx_cell++) {
//             for (int dy_cell = -1; dy_cell <= 1; dy_cell++) {

//                 int nx = (cx + dx_cell + nCells) % nCells;
//                 int ny = (cy + dy_cell + nCells) % nCells;

//                 int nc = cell_id(nx, ny);

//                 for (int j : cell_list[nc]) {

//                     if (j <= i) continue;

//                     // ========== YOUR ORIGINAL INTERACTION LOGIC ==========

//                     double dx = abs(P[j].x - P[i].x);
//                     double dy = abs(P[j].y - P[i].y);

//                     int dirx = 1, diry = 1;
//                     if(P[j].x > P[i].x){
//                         if(dx > L/2) dirx = 1;
//                         else dirx = -1;
//                     } else {
//                         if(dx > L/2) dirx = -1;
//                         else dirx = 1;
//                     }

//                     if(P[j].y > P[i].y){
//                         if(dy > L/2) diry = 1;
//                         else diry = -1;
//                     } else {
//                         if(dy > L/2) diry = -1;
//                         else diry = 1;
//                     }

//                     if(dx > L/2) dx = L - dx;
//                     if(dy > L/2) dy = L - dy;

//                     double d2 = dx*dx + dy*dy;
//                     if (d2 < 1e-12 || d2 >= r_int*r_int) continue;

//                     double d = sqrt(d2);
//                     double dxn = dx / d;
//                     double dyn = dy / d;

//                     double Fmag = beta_CF*(r_int/d - 1.0);
//                     double Fx = Fmag * dxn;
//                     double Fy = Fmag * dyn;

//                     Fx *= dirx;
//                     Fy *= diry;

//                     forces[i].first  += Fx;
//                     forces[i].second += Fy;
//                     forces[j].first  -= Fx;
//                     forces[j].second -= Fy;

//                     double qjx = cos(P[j].theta_n);
//                     double qjy = sin(P[j].theta_n);

//                     double weight = 0.5*(1.0 + dxn*qjx + dyn*qjy);

//                     JCF[i].first  += weight * dxn;
//                     JCF[i].second += weight * dyn;

//                     // ===================================================
//                 }
//             }}
//         }
//     }}

//     // ---------------- velocity + force output ----------------
//     for (int i = 0; i < N; i++) {

//         forceout << step << ","
//                  << P[i].x << ","
//                  << P[i].y << ","
//                  << P[i].vx << ","
//                  << P[i].vy << ","
//                  << P[i].theta_n << ","
//                  << forces[i].first << ","
//                  << forces[i].second << "\n";

//         P[i].vx = v0*cos(P[i].theta_n) + forces[i].first;
//         P[i].vy = v0*sin(P[i].theta_n) + forces[i].second;

//         P[i].fx = forces[i].first;
//         P[i].fy = forces[i].second;
//     }

//     // ---------------- position update ----------------
//     for (int i = 0; i < N; i++) {
//         P[i].x += P[i].vx * dt;
//         P[i].y += P[i].vy * dt;
//         wrap_position(P[i].x, P[i].y);
//     }

//     // ---------------- angle update ----------------
//     for (int i = 0; i < N; i++) {

//         double nx = cos(P[i].theta_n);
//         double ny = sin(P[i].theta_n);

//         double npx = -ny;
//         double npy =  nx;

//         P[i].theta_n += dt * alpha_CF *
//             (JCF[i].first*npx + JCF[i].second*npy);
//     }
// }



// ================= MAIN ======================
int main() {

    mt19937 rng(random_device{}());

    // box size from density
    L = sqrt(N / rho);

    cell_size = r_int;
    nCells = int(L / cell_size);
    if (nCells < 1) nCells = 1;
    //vector<vector<int>> cell_list(nCells * nCells);


    vector<Particle> P = initialize(N, rng);            

    ofstream fout("simulation_data.txt");
    fout << "# L=" << L << "\n";
    

    forceout.open("particle1_full_trace.txt");
    forceout << "# step x y vx vy theta Fx Fy\n";

    int frame=0;
    for (int step=0; step<STEPS; step++) {
        
        update(P, rng, step);

        if (step % SAVE_FREQ == 0)
            save_snapshot(fout, P, frame++);
    }

    fout.close();
    forceout.close();
    cout << "Simulation complete.\n";
    return 0;
}
