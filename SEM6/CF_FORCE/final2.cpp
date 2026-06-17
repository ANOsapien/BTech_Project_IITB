#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>

#include <cmath>
using namespace std;
ofstream forceout;
// ================= PARAMETERS =================
const int    N =100;

const double v0 =1;

const double r_int = 1.0;

const double beta_CF  = 0.0;
const double alpha_CF = 1.0;

const double Dr = 0.0;

const double rho = 0.01;

const double dt = 0.001;
const int    STEPS = 80000;
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

    // columns: x y vx vy theta fx fy
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

//lattice
// vector<Particle> initialize(int N, mt19937 &rng) {

//     vector<Particle> P(N);

//     int nSide = ceil(sqrt(N));      // lattice points per side
//     double a = L / nSide;           // lattice spacing

//     uniform_real_distribution<double> ang(0.0, 2*M_PI);

//     int idx = 0;
//     for (int i = 0; i < nSide && idx < N; i++) {
//         for (int j = 0; j < nSide && idx < N; j++) {

//             P[idx].x = (i + 0.5) * a;
//             P[idx].y = (j + 0.5) * a;

//             wrap_position(P[idx].x, P[idx].y);

//             P[idx].theta_n = ang(rng);
//             P[idx].vx = v0 * cos(P[idx].theta_n);
//             P[idx].vy = v0 * sin(P[idx].theta_n);

//             idx++;
//         }
//     }

//     return P;
// }


//circle 
// vector<Particle> initialize(int N, mt19937 &rng) {

//     vector<Particle> P(N);

//     double R  = 0.35 * L;     // ring radius (paper uses ~0.3–0.4 of box)
//     double cx = 0.5 * L;
//     double cy = 0.5 * L;

//     for (int i = 0; i < N; i++) {

//         double phi = 2.0 * M_PI * i / N;

//         // positions on ring
//         P[i].x = cx + R * cos(phi);
//         P[i].y = cy + R * sin(phi);
//         wrap_position(P[i].x, P[i].y);

//         // polarity tangential to ring (spiral/ring motion)
//         P[i].theta_n = phi + M_PI/2.0;

//         // overdamped velocity will be set in update(), but initialize anyway
//         P[i].vx = v0 * cos(P[i].theta_n);
//         P[i].vy = v0 * sin(P[i].theta_n);
//     }

//     return P;
// }

//circle touching always
// vector<Particle> initialize(int N, mt19937 &rng) {

//     vector<Particle> P(N);

//     double R = r_int / (2.0 * sin(M_PI / N))*0.5;
//     double cx = 0.5 * L;
//     double cy = 0.5 * L;

//     for (int i = 0; i < N; i++) {

//         double phi = 2.0 * M_PI * i / N +M_PI;

//         // positions on ring
//         P[i].x = cx + R * cos(phi);
//         P[i].y = cy + R * sin(phi);
//         wrap_position(P[i].x, P[i].y);

//         // polarity tangential to ring (spiral/ring motion)
//         P[i].theta_n = phi + M_PI/2 ;

//         // overdamped velocity will be set in update(), but initialize anyway
//         P[i].vx = v0 * cos(P[i].theta_n);
//         P[i].vy = v0 * sin(P[i].theta_n);
//     }

//     return P;
// }



//concentric rings
// vector<Particle> initialize(int N, mt19937 &rng){

//     vector<Particle> P;
//     P.reserve(N);

//     double cx = 0.5 * L;
//     double cy = 0.5 * L;

//     double compression = 0.5;   // <1 gives overlap

//     int count = 0;
//     int k = 1;

//     while(count < N){

//         double R = k * r_int * compression;

//         int Nk = int(round(2.0 * M_PI * k));
//         if(Nk < 6) Nk = 6;

//         for(int i=0;i<Nk && count < N;i++){

//             double phi = 2.0 * M_PI * i / Nk;

//             Particle p;

//             p.x = cx + R * cos(phi);
//             p.y = cy + R * sin(phi);
//             wrap_position(p.x, p.y);

//             p.theta_n = phi + M_PI/2.0;
//             p.vx = v0 * cos(p.theta_n);
//             p.vy = v0 * sin(p.theta_n);

//             P.push_back(p);
//             count++;
//         }
//         k++;
//     }

//     return P;
// }

//line
// vector<Particle> initialize(int N, mt19937 &rng) {

//     vector<Particle> P(N);

//     double cy = 0.5 * L;

//     double alpha = 0.45;              // <--- compression factor
//     double spacing = alpha * 2.0 * r_int;

//     double occupied = spacing * (N - 1);
//     double x0 = 0.5 * (L - occupied);

//     double theta = M_PI/2;

//     for (int i = 0; i < N; i++) {

//         P[i].x = x0 + i * spacing;
//         P[i].y = cy;

//         wrap_position(P[i].x, P[i].y);

//         P[i].theta_n = theta;

//         P[i].vx = v0 * cos(theta);
//         P[i].vy = v0 * sin(theta);
//     }

//     return P;
// }


//line with angle bias
// vector<Particle> initialize(int N, mt19937 &rng) {

//     vector<Particle> P(N);

//     double cy = 0.5 * L;

//     double alpha = 0.2;
//     double spacing = alpha * 2.0 * r_int;

//     double occupied = spacing * (N - 1);
//     double x0 = 0.5 * (L - occupied);

//     for (int i = 0; i < N; i++) {

//         P[i].x = x0 + i * spacing;
//         P[i].y = cy;

//         wrap_position(P[i].x, P[i].y);

//         if (i == 0)
//             P[i].theta_n =  0.79;     // +90°
//         else if (i == N-1)
//             P[i].theta_n = -3.14;     // -90°
//         else
//             P[i].theta_n = 0.0;         // straight

//         P[i].vx = v0 * cos(P[i].theta_n);
//         P[i].vy = v0 * sin(P[i].theta_n);
//     }

//     return P;
// }


//2 particle bias
// vector<Particle> initialize(int N, mt19937 &rng) {

//     vector<Particle> P(N);

//     double cy = 0.5 * L;

//     double alpha = 0.2;
//     double spacing = alpha * 2.0 * r_int;

//     double occupied = spacing * (N - 1);
//     double x0 = 0.5 * (L - occupied);

//     for (int i = 0; i < N; i++) {

//         P[i].x = x0 + i * spacing;
//         P[i].y = cy;

//         wrap_position(P[i].x, P[i].y);

//         if (N == 2) {
//             if (i == 0)
//                 P[i].theta_n = 0;   // 45°
//             else
//                 P[i].theta_n = 0;        // 0°
//         }
//         else {
//             P[i].theta_n = 0.0;            // default for larger N
//         }

//         P[i].vx = v0 * cos(P[i].theta_n);
//         P[i].vy = v0 * sin(P[i].theta_n);
//     }

//     return P;
// }


// vector<Particle> initialize(int N, mt19937 &rng) {

//     vector<Particle> P(N);

//     double R0 = 0.35 * L;      // radius of dense cluster
//     double cx = 0.5 * L;
//     double cy = 0.5 * L;

//     uniform_real_distribution<double> uni(0.0,1.0);
//     uniform_real_distribution<double> ang(0.0,2*M_PI);

//     for(int i=0;i<N;i++){

//         // uniform filling of a disk (not ring!)
//         double r = R0 * sqrt(uni(rng));
//         double phi = ang(rng);

//         P[i].x = cx + r*cos(phi);
//         P[i].y = cy + r*sin(phi);
//         wrap_position(P[i].x, P[i].y);

//         // polarity initially random (paper does this)
//         P[i].theta_n = ang(rng);

//         P[i].vx = v0*cos(P[i].theta_n);
//         P[i].vy = v0*sin(P[i].theta_n);
//     }

//     return P;
// }

// vector<Particle> initialize(int N) {

//     vector<Particle> P(N);

//     double R = 0.35 * L;          // ring radius
//     double cx = 0.5 * L;
//     double cy = 0.5 * L;

//     for (int i = 0; i < N; i++) {

//         double phi = 2.0 * M_PI * i / N;

//         P[i].x = cx + R * cos(phi);
//         P[i].y = cy + R * sin(phi);
//         wrap_position(P[i].x, P[i].y);

//         // tangential direction
//         P[i].theta_n = phi + M_PI/4;

//         P[i].vx = -v0 * cos(P[i].theta_n);
//         P[i].vy = -v0 * sin(P[i].theta_n);
//     }

//     return P;
// }
//circle with distance bias
vector<Particle> initialize(int N, mt19937 &rng) {

    vector<Particle> P(N);

    double cx = 0.5 * L;
    double cy = 0.5 * L;

    // perfect touching ring radius
    double R = r_int / (2.0 * sin(M_PI / N));

    uniform_real_distribution<double> noise(-0.01*r_int, 0.01*r_int);

    for(int i=0;i<N;i++){

        double phi = 2.0*M_PI*i/N ;

        // perfect circle
        double x0 = cx + R*cos(phi);
        double y0 = cy + R*sin(phi);

        // small distortion (keeps near touching)
        P[i].x = x0 + noise(rng);
        P[i].y = y0 + noise(rng);

        wrap_position(P[i].x,P[i].y);

        // tangential polarity
        P[i].theta_n = phi + M_PI/2;

        P[i].vx = v0*cos(P[i].theta_n);
        P[i].vy = v0*sin(P[i].theta_n);
    }

    return P;
}


bool nearly_equal(double a, double b,
                  double abs_tol = 1e-12,
                  double rel_tol = 1e-9)
{
    return std::fabs(a - b) <= std::max(abs_tol, rel_tol * std::max(std::fabs(a), std::fabs(b)));
}


// void update(vector<Particle>& P, mt19937 &rng, int step){
//     int N = P.size();


//     // FIX 1: Initialize forces and JCF to {0,0} to avoid garbage values
//     vector<pair<double,double>> forces(N, {0.0, 0.0});
//     vector<pair<double,double>> JCF(N,    {0.0, 0.0});

//     for(int i=0; i<N; i++){
//         for(int j=i+1; j<N; j++){

//             double dx = P[j].x - P[i].x;   // signed
//             double dy = P[j].y - P[i].y;   // signed

//             // PBC: wrap to [-L/2, L/2]
//             if (dx >  L/2) dx -= L;
//             if (dx < -L/2) dx += L;
//             if (dy >  L/2) dy -= L;
//             if (dy < -L/2) dy += L;

//             double d2 = dx*dx + dy*dy;
//             if (d2 < 1e-12 || d2 >= r_int*r_int) continue;

//             double d   = sqrt(d2);
//             double dxn = dx / d;   // unit vector from i toward j
//             double dyn = dy / d;

//             // Repulsion: push i away from j, push j away from i
//             double Fmag = beta_CF * (r_int/d - 1.0);
//             forces[i].first  -= Fmag * dxn;
//             forces[i].second -= Fmag * dyn;
//             forces[j].first  += Fmag * dxn;
//             forces[j].second += Fmag * dyn;

//             // FIX 2: JCF weight: was dyn*dyn*qjy (typo), now dyn*qjy (correct dot product)
//             double qjx   = cos(P[j].theta_n);
//             double qjy   = sin(P[j].theta_n);
//             double weight = 0.5*(1.0 + dxn*qjx + dyn*qjy);  // FIXED: was dyn*dyn*qjy
//             JCF[i].first  += weight * dxn;
//             JCF[i].second += weight * dyn;
//         }
//     }

//     // FIX 3: Use double (not float) for total force accumulation
//     double tfx = 0.0, tfy = 0.0;

//     for (int i = N-1; i >=0; i--) {
//         const Particle& p = P[i];

//         // FIX 4: Log forces[i] (not always forces[0])
//         forceout << step << ","
//                 << p.x << ","
//                 << p.y << ","
//                 << p.vx << ","
//                 << p.vy << ","
//                 << p.theta_n << ","
//                 << forces[i].first << ","    // FIXED: was forces[0].first
//                 << forces[i].second << "\n"; // FIXED: was forces[0].second

//         P[i].vx = v0 * cos(P[i].theta_n) + forces[i].first;
//         P[i].vy = v0 * sin(P[i].theta_n) + forces[i].second;
//         P[i].fx = forces[i].first;
//         P[i].fy = forces[i].second;

//         tfx += P[i].fx;
//         tfy += P[i].fy;
//     }

//     // ---- position update ----
//     for (int i = N-1; i >=0; i--) {
//         P[i].x += P[i].vx * dt;
//         P[i].y += P[i].vy * dt;
//         wrap_position(P[i].x, P[i].y);
//     }

//     // ---- angle update ----
//     for (int i = N-1; i >=0; i--) {
//         double nx = cos(P[i].theta_n);
//         double ny = sin(P[i].theta_n);
//         double npx = -ny;
//         double npy =  nx;

//         P[i].theta_n += dt * alpha_CF *
//             (JCF[i].first*npx + JCF[i].second*npy);
//     }

//     // FIX 5: Condition was inverted — total force SHOULD be ~zero (Newton's 3rd law)
//     // Warn only when it is NOT zero
//     if(!nearly_equal(tfx, 0.0) || !nearly_equal(tfy, 0.0)){
//         cout << "Force imbalance at step " << step
//              << "  tfx=" << tfx << "  tfy=" << tfy << endl;
//     }
// }


void update(vector<Particle>& P, mt19937 &rng, int step){

    int N = P.size();

    vector<pair<double,double>> forces(N,{0.0,0.0});
    vector<pair<double,double>> JCF(N,{0.0,0.0});

    // ---------- interaction forces & torques ----------
    for(int i=0;i<N;i++){
        for(int j=i+1;j<N;j++){

            double dx = P[j].x - P[i].x;
            double dy = P[j].y - P[i].y;

            if(dx >  L/2) dx -= L;
            if(dx < -L/2) dx += L;
            if(dy >  L/2) dy -= L;
            if(dy < -L/2) dy += L;

            double d2 = dx*dx + dy*dy;
            if(d2 < 1e-12 || d2 >= r_int*r_int) continue;

            double d = sqrt(d2);
            double dxn = dx/d;
            double dyn = dy/d;

            double Fmag = beta_CF * (r_int/d - 1.0);

            forces[i].first  -= Fmag*dxn;
            forces[i].second -= Fmag*dyn;
            forces[j].first  += Fmag*dxn;
            forces[j].second += Fmag*dyn;

            double qjx = cos(P[j].theta_n);
            double qjy = sin(P[j].theta_n);

            double w = 0.5*(1.0 + dxn*qjx + dyn*qjy);

            JCF[i].first  += w*dxn;
            JCF[i].second += w*dyn;
        }
    }

    // ---------- new state containers ----------
    vector<double> newx(N), newy(N), newtheta(N);
    vector<double> newvx(N), newvy(N);

    // ---------- integrate from OLD state ----------
    for(int i=0;i<N;i++){

        newvx[i] = v0*cos(P[i].theta_n) + forces[i].first;
        newvy[i] = v0*sin(P[i].theta_n) + forces[i].second;

        newx[i] = P[i].x + newvx[i]*dt;
        newy[i] = P[i].y + newvy[i]*dt;

        wrap_position(newx[i], newy[i]);

        double nx = cos(P[i].theta_n);
        double ny = sin(P[i].theta_n);

        double npx = -ny;
        double npy =  nx;

        newtheta[i] = P[i].theta_n
            + dt*alpha_CF*(JCF[i].first*npx + JCF[i].second*npy);
    }

    // ---------- commit simultaneously ----------
    for(int i=0;i<N;i++){
        P[i].x = newx[i];
        P[i].y = newy[i];
        P[i].theta_n = newtheta[i];
        P[i].vx = newvx[i];
        P[i].vy = newvy[i];
        P[i].fx = forces[i].first;
        P[i].fy = forces[i].second;
    }
}


// ================= MAIN ======================
int main() {

    mt19937 rng(random_device{}());

    // box size from density
    L = sqrt(N / rho);

    cell_size = r_int;
    nCells = int(L / cell_size);
    if (nCells < 1) nCells = 1;

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