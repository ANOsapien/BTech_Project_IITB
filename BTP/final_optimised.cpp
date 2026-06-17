#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>
#include <cmath>
using namespace std;

// ================= PARAMETERS =================
const int    N    = 2;
const double v0   = 0.1;
const double r_int = 1.0;

const double beta_CF  = 0;
const double alpha_CF = 1.0;

const double Dr  = 0;
const double rho = 0.01;

const double dt        = 0.001;
const int    STEPS     = 40000;
const int    SAVE_FREQ = 10;

// ================= NEIGHBOR LIST ==============
const double r_skin = 0.3;
const double r_cut  = r_int + r_skin;
const double r_cut2 = r_cut * r_cut;           // precomputed
const double r_int2 = r_int * r_int;           // precomputed
const double r_skin2_quarter = 0.25 * r_skin * r_skin; // precomputed

vector<vector<int>> neighbors;
vector<double> lastx, lasty;

// ================= GLOBALS ====================
double L;
double halfL;   // precomputed L/2

// ================= PARTICLE ==================
struct Particle {
    double x, y;
    double theta_n;
    double cos_theta, sin_theta;   // cached trig — updated whenever theta changes
    double vx, vy;
    double fx, fy;
};

// ================= PBC =======================
inline void wrap_position(double &x, double &y){
    if(x <  0) x += L;
    if(x >= L) x -= L;
    if(y <  0) y += L;
    if(y >= L) y -= L;
}

inline void wrap_distance(double &dx, double &dy){
    if(dx >  halfL) dx -= L;
    if(dx < -halfL) dx += L;
    if(dy >  halfL) dy -= L;
    if(dy < -halfL) dy += L;
}

//===================ANGLE CAPPING===============
inline double wrap_angle(double theta){
    theta = fmod(theta, 2.0*M_PI);
    if(theta < 0) theta += 2.0*M_PI;
    return theta;
}

// ================= SNAPSHOT ==================
void save_snapshot(ofstream &fout, const vector<Particle>& P, int frame){
    fout << "FRAME " << frame << "\n";
    for(auto &p : P)
        fout << p.x  << " " << p.y  << " "
             << p.vx << " " << p.vy << " "
             << p.theta_n << " " << p.fx << " " << p.fy << "\n";
}

// ================= INITIALIZE =================
//2 particles
vector<Particle> initialize(int N, mt19937 &rng){

    vector<Particle> P(N);

    double cy      = 0.5 * L;
    double alpha   = 0.4;
    double spacing = alpha * 2.0 * r_int;
    double x0      = 0.5 * (L - spacing * (N - 1));

    for(int i = 0; i < N; i++){
        P[i].x = x0 + i * spacing;
        P[i].y = cy;
        wrap_position(P[i].x, P[i].y);

        if(i == 0)      P[i].theta_n = 0;
        else            P[i].theta_n = 0;

        P[i].cos_theta = cos(P[i].theta_n);
        P[i].sin_theta = sin(P[i].theta_n);

        P[i].vx = v0 * P[i].cos_theta;
        P[i].vy = v0 * P[i].sin_theta;
        P[i].fx = P[i].fy = 0.0;
    }
    return P;
}

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


// concentric ring with velocity
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

//         for(int i = 0; i < Nk && count < N; i++){

//             double phi = 2.0 * M_PI * i / Nk;

//             Particle p;

//             p.x = cx + R * cos(phi);
//             p.y = cy + R * sin(phi);
//             wrap_position(p.x, p.y);

//             // tangential direction = phi + pi/2 (CCW)
//             p.theta_n   = phi + M_PI / 2.0;
//             p.cos_theta = cos(p.theta_n);
//             p.sin_theta = sin(p.theta_n);

//             // velocity aligned with theta_n — tangential to ring
//             p.vx = v0 * p.cos_theta;   // = -v0 * sin(phi)
//             p.vy = v0 * p.sin_theta;   // =  v0 * cos(phi)

//             p.fx = p.fy = 0.0;

//             P.push_back(p);
//             count++;
//         }
//         k++;
//     }

//     return P;
// }


//out velocity
// vector<Particle> initialize(int N, mt19937 &rng) {

//     vector<Particle> P;
//     P.reserve(N);

//     const double cx = 0.5 * L;
//     const double cy = 0.5 * L;

//     const double compression = 0.5;
//     const double omega = 1.0;

//     int count = 0;
//     int k = 1;

//     while (count < N) {

//         double R = k * r_int * compression;

//         int Nk = int(round(2.0 * M_PI * k));
//         if (Nk < 6) Nk = 6;

//         int ring_start = count;

//         for (int i = 0; i < Nk && count < N; i++) {

//             double phi = 2.0 * M_PI * i / Nk;

//             Particle p;

//             p.x = cx + R * cos(phi);
//             p.y = cy + R * sin(phi);

//             p.theta_n = phi + M_PI / 2.0;

//             // No velocity initially
//             p.vx = 0.0;
//             p.vy = 0.0;

//             wrap_position(p.x, p.y);

//             P.push_back(p);
//             count++;
//         }

//         // If this ring completed the particle set → outermost ring
//         if (count == N) {

//             for (int j = ring_start; j < count; j++) {

//                 double dx = P[j].x - cx;
//                 double dy = P[j].y - cy;

//                 double phi = atan2(dy, dx);
//                 double R_ring = sqrt(dx*dx + dy*dy);

//                 double theta = phi + M_PI / 2.0;

//                 double v = omega * R_ring;

//                 P[j].vx = v * cos(theta);
//                 P[j].vy = v * sin(theta);
//             }
//         }

//         k++;
//     }

//     return P;
// }



//in velocity
// vector<Particle> initialize(int N, mt19937 &rng) {

//     vector<Particle> P;
//     P.reserve(N);

//     const double cx = 0.5 * L;
//     const double cy = 0.5 * L;

//     const double compression = 1.1;
//     const double omega = 1.0;

//     int count = 0;
//     int k = 1;

//     while (count < N) {

//         double R = k * r_int * compression;

//         int Nk = int(round(2.0 * M_PI * k));
//         if (Nk < 6) Nk = 6;

//         for (int i = 0; i < Nk && count < N; i++) {

//             double phi = 2.0 * M_PI * i / Nk;

//             Particle p;

//             p.x = cx + R * cos(phi);
//             p.y = cy + R * sin(phi);

//             p.theta_n = phi + M_PI / 2.0;

//             if (k == 1) {   // only innermost ring moves
//                 double v = omega * R;
//                 p.vx = v * cos(p.theta_n);
//                 p.vy = v * sin(p.theta_n);
//             } 
//             else {
//                 p.vx = 0.0;
//                 p.vy = 0.0;
//             }

//             wrap_position(p.x, p.y);

//             P.push_back(p);
//             count++;
//         }

//         k++;
//     }

//     return P;
// }


// half circular bunch moving from each side
// vector<Particle> initialize(int N, mt19937 &rng){

//     vector<Particle> P;
//     P.reserve(N);

//     int N_half = N / 2;

//     double r_cluster = 5.0 * r_int;   // radius of each bunch
//     double cx_left  = 0.2 * L;
//     double cx_right = 0.8 * L;
//     double cy = 0.5 * L;

//     uniform_real_distribution<double> ang_dist(0.0, 2*M_PI);
//     uniform_real_distribution<double> rad_dist(0.0, r_cluster);

//     // Left cluster (moving right)
//     for(int i=0;i<N_half;i++){

//         double a = ang_dist(rng);
//         double r = sqrt(rad_dist(rng)*rad_dist(rng));

//         Particle p;

//         p.x = cx_left + r*cos(a);
//         p.y = cy + r*sin(a);
//         wrap_position(p.x,p.y);

//         p.theta_n = 0.0;     // direction → right
//         p.vx = v0;
//         p.vy = 0.0;

//         P.push_back(p);
//     }

//     // Right cluster (moving left)
//     for(int i=N_half;i<N;i++){

//         double a = ang_dist(rng);
//         double r = sqrt(rad_dist(rng)*rad_dist(rng));

//         Particle p;

//         p.x = cx_right + r*cos(a);
//         p.y = cy + r*sin(a);
//         wrap_position(p.x,p.y);

//         p.theta_n = M_PI;    // direction → left
//         p.vx = -v0;
//         p.vy = 0.0;

//         P.push_back(p);
//     }

//     return P;
// }

// vertical lines moving
// vector<Particle> initialize(int N, mt19937 &rng){

//     vector<Particle> P;
//     P.reserve(N);

//     int N_half = N/2;

//     double x_left  = 0.2 * L;
//     double x_right = 0.8 * L;

//     double y_spacing = L / (N_half + 1);

//     // Left line → moving right
//     for(int i=0;i<N_half;i++){

//         Particle p;

//         p.x = x_left;
//         p.y = (i+1) * y_spacing;
//         wrap_position(p.x,p.y);

//         p.theta_n = 0.0;     // toward right
//         p.vx = v0;
//         p.vy = 0.0;

//         P.push_back(p);
//     }

//     // Right line → moving left
//     for(int i=N_half;i<N;i++){

//         Particle p;

//         p.x = x_right;
//         p.y = (i - N_half + 1) * y_spacing;
//         wrap_position(p.x,p.y);

//         p.theta_n = M_PI;    // toward left
//         p.vx = -v0;
//         p.vy = 0.0;

//         P.push_back(p);
//     }

//     return P;
// }


//horizontal rows moving 
// vector<Particle> initialize(int N, mt19937 &rng){

//     vector<Particle> P;
//     P.reserve(N);

//     int per_side = N/2;

//     int slabs = 6;                 // number of horizontal bands
//     int rows_per_slab = 6;         // thickness of each band

//     double dx = 2.1 * r_int;       // particle spacing (almost touching)
//     double dy = 2.1 * r_int;

//     double slab_gap = 4 * r_int;

//     double x_left_start  = 0.05 * L;
//     double x_right_start = 0.95 * L;

//     int count = 0;

//     for(int s=0; s<slabs && count < per_side; s++){

//         double y0 = (s+1)*slab_gap + s*rows_per_slab*dy;

//         for(int r=0; r<rows_per_slab && count < per_side; r++){

//             double y = y0 + r*dy;

//             for(int i=0; count < per_side; i++){

//                 double x = x_left_start + i*dx;
//                 if(x > 0.4*L) break;

//                 Particle p;

//                 p.x = x;
//                 p.y = y;
//                 wrap_position(p.x,p.y);

//                 p.theta_n = 0.0;
//                 p.vx = v0;
//                 p.vy = 0.0;

//                 P.push_back(p);
//                 count++;
//             }
//         }
//     }

//     for(int s=0; s<slabs && count < N; s++){

//         double y0 = (s+1)*slab_gap + s*rows_per_slab*dy;

//         for(int r=0; r<rows_per_slab && count < N; r++){

//             double y = y0 + r*dy;

//             for(int i=0; count < N; i++){

//                 double x = x_right_start - i*dx;
//                 if(x < 0.6*L) break;

//                 Particle p;

//                 p.x = x;
//                 p.y = y;
//                 wrap_position(p.x,p.y);

//                 p.theta_n = M_PI;
//                 p.vx = -v0;
//                 p.vy = 0.0;

//                 P.push_back(p);
//                 count++;
//             }
//         }
//     }

//     return P;
// }



//two rows 
// vector<Particle> initialize(int N, mt19937 &rng) {

//     vector<Particle> P;
//     P.reserve(N);

//     int N_half = N/2;

//     int nSide = ceil(sqrt(N_half));
//     double gap = 0.2 * L;

//     double width_block = (L - gap) / 2.0;
//     double a = width_block / nSide;

//     double x_left_start  = 0.0;
//     double x_right_start = width_block + gap;

//     int idx = 0;

//     // LEFT lattice (velocity → right)
//     for(int i=0;i<nSide && idx<N_half;i++){
//         for(int j=0;j<nSide && idx<N_half;j++){

//             Particle p;

//             p.x = x_left_start + (i + 0.5)*a;
//             p.y = (j + 0.5)*a;

//             wrap_position(p.x,p.y);

//             p.theta_n = 0.0;
//             p.vx = v0;
//             p.vy = 0.0;

//             P.push_back(p);
//             idx++;
//         }
//     }

//     // RIGHT lattice (velocity → left)
//     for(int i=0;i<nSide && idx<N;i++){
//         for(int j=0;j<nSide && idx<N;j++){

//             Particle p;

//             p.x = x_right_start + (i + 0.5)*a;
//             p.y = (j + 0.5)*a;

//             wrap_position(p.x,p.y);

//             p.theta_n = M_PI;
//             p.vx = -v0;
//             p.vy = 0.0;

//             P.push_back(p);
//             idx++;
//         }
//     }

//     return P;
// }

// 4 lattice blocks
// vector<Particle> initialize(int N, mt19937 &rng) {

//     vector<Particle> P;
//     P.reserve(N);

//     int per_block = N / 4;

//     int nSide = ceil(sqrt(per_block));

//     double gap_x = 0.2 * L;
//     double gap_y = 0.2 * L;

//     double block_w = (L - gap_x) / 2.0;
//     double block_h = (L - gap_y) / 2.0;

//     double a = min(block_w, block_h) / nSide;

//     double x_left  = 0.0;
//     double x_right = block_w + gap_x;

//     double y_bottom = 0.0;
//     double y_top    = block_h + gap_y;

//     int idx = 0;

//     // bottom-left block (→)
//     for(int i=0;i<nSide && idx<per_block;i++){
//         for(int j=0;j<nSide && idx<per_block;j++){

//             Particle p;

//             p.x = x_left + (i+0.5)*a;
//             p.y = y_bottom + (j+0.5)*a;

//             wrap_position(p.x,p.y);

//             p.theta_n = 0.0;
//             p.vx = v0;
//             p.vy = 0.0;

//             P.push_back(p);
//             idx++;
//         }
//     }

//     // top-left block (→)
//     for(int i=0;i<nSide && idx<2*per_block;i++){
//         for(int j=0;j<nSide && idx<2*per_block;j++){

//             Particle p;

//             p.x = x_left + (i+0.5)*a;
//             p.y = y_top + (j+0.5)*a;

//             wrap_position(p.x,p.y);

//             p.theta_n = 0.0;
//             p.vx = v0;
//             p.vy = 0.0;

//             P.push_back(p);
//             idx++;
//         }
//     }

//     // bottom-right block (←)
//     for(int i=0;i<nSide && idx<3*per_block;i++){
//         for(int j=0;j<nSide && idx<3*per_block;j++){

//             Particle p;

//             p.x = x_right + (i+0.5)*a;
//             p.y = y_bottom + (j+0.5)*a;

//             wrap_position(p.x,p.y);

//             p.theta_n = M_PI;
//             p.vx = -v0;
//             p.vy = 0.0;

//             P.push_back(p);
//             idx++;
//         }
//     }

//     // top-right block (←)
//     for(int i=0;i<nSide && idx<N;i++){
//         for(int j=0;j<nSide && idx<N;j++){

//             Particle p;

//             p.x = x_right + (i+0.5)*a;
//             p.y = y_top + (j+0.5)*a;

//             wrap_position(p.x,p.y);

//             p.theta_n = M_PI;
//             p.vx = -v0;
//             p.vy = 0.0;

//             P.push_back(p);
//             idx++;
//         }
//     }

//     return P;
// }

// trinagular lattice
// vector<Particle> initialize(int N, mt19937 &rng){

//     vector<Particle> P;
//     P.reserve(N);

//     int sites = N / 3;                    // 3 particles per triangle

//     double a = L / sqrt(sites);           // approximate lattice spacing
//     double dy = a * sqrt(3.0)/2.0;        // hex vertical spacing

//     double r_tri = 0.4 * r_int;           // triangle size

//     uniform_real_distribution<double> ang(0.0, 2*M_PI);

//     int idx = 0;
//     int row = 0;

//     for(double y = dy/2; y < L && idx < N; y += dy, row++){

//         double x_offset = (row % 2) ? a/2 : 0;

//         for(double x = a/2 + x_offset; x < L && idx < N; x += a){

//             // velocity chosen per triangle (lattice site)
//             double theta_site = ang(rng);
//             double vx_site = v0 * cos(theta_site);
//             double vy_site = v0 * sin(theta_site);

//             for(int k=0;k<3 && idx<N;k++){

//                 double phi = 2*M_PI*k/3.0;

//                 Particle p;

//                 p.x = x + r_tri*cos(phi);
//                 p.y = y + r_tri*sin(phi);

//                 wrap_position(p.x,p.y);

//                 p.theta_n = theta_site;
//                 p.vx = vx_site;
//                 p.vy = vy_site;

//                 P.push_back(p);
//                 idx++;
//             }
//         }
//     }

//     return P;
// }

// hexagonal lattice
// vector<Particle> initialize(int N, mt19937 &rng){

//     vector<Particle> P;
//     P.reserve(N);

//     int sites = N / 6;                 // 6 particles per hexagon

//     double a = L / sqrt(sites);        // approximate lattice spacing
//     double dy = a * sqrt(3.0)/2.0;     // hex lattice vertical spacing

//     double r_hex = 0.4 * r_int;        // hexagon size

//     uniform_real_distribution<double> ang(0.0, 2*M_PI);

//     int idx = 0;
//     int row = 0;

//     for(double y = dy/2; y < L && idx < N; y += dy, row++){

//         double x_offset = (row % 2) ? a/2 : 0;

//         for(double x = a/2 + x_offset; x < L && idx < N; x += a){

//             // velocity for this hexagon site
//             double theta_site = ang(rng);
//             double vx_site = v0 * cos(theta_site);
//             double vy_site = v0 * sin(theta_site);

//             for(int k=0;k<6 && idx<N;k++){

//                 double phi = 2*M_PI*k/6.0;

//                 Particle p;

//                 p.x = x + r_hex*cos(phi);
//                 p.y = y + r_hex*sin(phi);

//                 wrap_position(p.x,p.y);

//                 p.theta_n = theta_site;
//                 p.vx = vx_site;
//                 p.vy = vy_site;

//                 P.push_back(p);
//                 idx++;
//             }
//         }
//     }

//     return P;
// }

// shear flow
// vector<Particle> initialize(int N, mt19937 &rng){

//     vector<Particle> P(N);

//     int nSide = ceil(sqrt(N));
//     double a = L / nSide;        // lattice spacing

//     double gamma = 1.0;          // shear rate

//     int idx = 0;

//     for(int i=0;i<nSide && idx<N;i++){
//         for(int j=0;j<nSide && idx<N;j++){

//             double x = (i + 0.5)*a;
//             double y = (j + 0.5)*a;

//             P[idx].x = x;
//             P[idx].y = y;

//             wrap_position(P[idx].x, P[idx].y);

//             double vx = gamma * (y - 0.5*L);
//             double vy = 0.0;

//             P[idx].vx = vx;
//             P[idx].vy = vy;

//             P[idx].theta_n = atan2(vy, vx);

//             idx++;
//         }
//     }

//     return P;
// }


//4 lines
// vector<Particle> initialize(int N, mt19937 &rng) {

//     vector<Particle> P;
//     P.reserve(N);

//     // ── TUNE THESE ──────────────────────────────
//     const double spacing_x = 0.5 * r_int;  // horizontal compression
//     const double spacing_y = 1.5 * r_int;  // vertical compression
//     const double offset_x  = 0.25 * L;     // how far blocks start from center
//     const double offset_y  = 0.15 * L;     // vertical separation between pairs
//     // ────────────────────────────────────────────

//     int per_block = N / 4;
//     int nSide     = (int)ceil(sqrt((double)per_block));

//     double cx = 0.5 * L;
//     double cy = 0.5 * L;

//     double bx[4]     = { cx - offset_x, cx + offset_x,
//                          cx - offset_x, cx + offset_x };

//     double by[4]     = { cy + offset_y, cy + offset_y,
//                          cy - offset_y, cy - offset_y };

//     double bvx[4]    = {  v0, -v0,  v0, -v0 };
//     double bvy[4]    = { 0.0, 0.0, 0.0, 0.0 };

//     double btheta[4] = { 0.0, M_PI, 0.0, M_PI };

//     for (int b = 0; b < 4; b++) {

//         int count = 0;
//         double x0 = bx[b] - 0.5 * (nSide - 1) * spacing_x;
//         double y0 = by[b] - 0.5 * (nSide - 1) * spacing_y;

//         for (int i = 0; i < nSide && count < per_block; i++) {
//             for (int j = 0; j < nSide && count < per_block; j++) {

//                 Particle p;

//                 p.x = x0 + i * spacing_x;
//                 p.y = y0 + j * spacing_y;

//                 wrap_position(p.x, p.y);

//                 p.theta_n   = btheta[b];
//                 p.cos_theta = cos(p.theta_n);
//                 p.sin_theta = sin(p.theta_n);
//                 p.vx        = bvx[b];
//                 p.vy        = bvy[b];
//                 p.fx = p.fy = 0.0;

//                 P.push_back(p);
//                 count++;
//             }
//         }
//     }

//     return P;
// }


// ================= NEIGHBOR LIST =================
void build_neighbor_list(const vector<Particle>& P){

    int n = P.size();
    neighbors.assign(n, {});

    for(int i = 0; i < n; i++)
        for(int j = i + 1; j < n; j++){
            double dx = P[j].x - P[i].x;
            double dy = P[j].y - P[i].y;
            wrap_distance(dx, dy);

            if(dx*dx + dy*dy < r_cut2){
                neighbors[i].push_back(j);
                neighbors[j].push_back(i);
            }
        }

    lastx.resize(n);
    lasty.resize(n);
    for(int i = 0; i < n; i++){
        lastx[i] = P[i].x;
        lasty[i] = P[i].y;
    }
}

bool need_rebuild(const vector<Particle>& P){
    for(int i = 0; i < (int)P.size(); i++){
        double dx = P[i].x - lastx[i];
        double dy = P[i].y - lasty[i];
        wrap_distance(dx, dy);
        if(dx*dx + dy*dy > r_skin2_quarter)
            return true;
    }
    return false;
}

// ================= FORCE COMPUTE =================
void compute_forces(const vector<Particle>& P,
                    vector<pair<double,double>>& forces,
                    vector<pair<double,double>>& JCF){

    int n = P.size();

    for(int i = 0; i < n; i++){
        forces[i] = {0.0, 0.0};
        JCF[i]    = {0.0, 0.0};
    }

    for(int i = 0; i < n; i++){
        for(int j : neighbors[i]){

            if(j <= i) continue;

            double dx = P[j].x - P[i].x;
            double dy = P[j].y - P[i].y;
            wrap_distance(dx, dy);

            double d2 = dx*dx + dy*dy;
            if(d2 < 1e-12 || d2 >= r_int2) continue;

            double inv_d = 1.0 / sqrt(d2);       // one sqrt, use reciprocal
            double dxn   = dx * inv_d;
            double dyn   = dy * inv_d;

            // conservative force
            if(beta_CF != 0.0){
                double Fmag = beta_CF * (r_int * inv_d - 1.0);
                forces[i].first  -= Fmag * dxn;
                forces[i].second -= Fmag * dyn;
                forces[j].first  += Fmag * dxn;
                forces[j].second += Fmag * dyn;
            }

            // non-reciprocal torque — use cached trig
            if(alpha_CF != 0.0){
                // torque on i from j
                double dot_j = dxn * P[j].cos_theta + dyn * P[j].sin_theta;
                double wi    = 0.5 * (1.0 + dot_j);
                JCF[i].first  += wi * dxn;
                JCF[i].second += wi * dyn;

                // torque on j from i
                double dot_i = dxn * P[i].cos_theta + dyn * P[i].sin_theta;
                double wj    = 0.5 * (1.0 - dot_i);
                JCF[j].first  -= wj * dxn;
                JCF[j].second -= wj * dyn;
            }
        }
    }
}

// ================= HEUN UPDATE =================
void update(vector<Particle>& P, mt19937 &rng){

    int n = P.size();

    // reuse static buffers to avoid repeated allocation
    static vector<pair<double,double>> F0, J0, F1, J1;
    static vector<Particle> Ppred;

    F0.resize(n); J0.resize(n);
    F1.resize(n); J1.resize(n);
    Ppred.resize(n);

    static normal_distribution<double> gauss(0.0,1.0);
    double noise_amp = sqrt(2.0 * Dr * dt);

    compute_forces(P, F0, J0);

    Ppred = P;

    // --- predictor ---
    for(int i = 0; i < n; i++){
        double vx0 = v0 * P[i].cos_theta + F0[i].first;
        double vy0 = v0 * P[i].sin_theta + F0[i].second;

        Ppred[i].x += vx0 * dt;
        Ppred[i].y += vy0 * dt;
        wrap_position(Ppred[i].x, Ppred[i].y);

        double dtheta = dt * alpha_CF * (J0[i].first  * (-P[i].sin_theta)
                                       + J0[i].second *   P[i].cos_theta);
        Ppred[i].theta_n += dtheta;

        // update cached trig for predicted state
        Ppred[i].cos_theta = cos(Ppred[i].theta_n);
        Ppred[i].sin_theta = sin(Ppred[i].theta_n);
    }

    compute_forces(Ppred, F1, J1);

    // --- corrector ---
    for(int i = 0; i < n; i++){
        double vx0 = v0 * P[i].cos_theta    + F0[i].first;
        double vy0 = v0 * P[i].sin_theta    + F0[i].second;
        double vx1 = v0 * Ppred[i].cos_theta + F1[i].first;
        double vy1 = v0 * Ppred[i].sin_theta + F1[i].second;

        P[i].x += 0.5 * (vx0 + vx1) * dt;
        P[i].y += 0.5 * (vy0 + vy1) * dt;
        wrap_position(P[i].x, P[i].y);

        double d0 = alpha_CF * (J0[i].first  * (-P[i].sin_theta)
                              + J0[i].second *   P[i].cos_theta);
        double d1 = alpha_CF * (J1[i].first  * (-Ppred[i].sin_theta)
                              + J1[i].second *   Ppred[i].cos_theta);


        double noise = noise_amp * gauss(rng);
        P[i].theta_n += 0.5 * (d0 + d1) * dt + noise;

        P[i].theta_n = wrap_angle(P[i].theta_n);


        // update cached trig for corrected state
        P[i].cos_theta = cos(P[i].theta_n);
        P[i].sin_theta = sin(P[i].theta_n);

        P[i].vx = vx1;
        P[i].vy = vy1;
        P[i].fx = F1[i].first;
        P[i].fy = F1[i].second;
    }
}

// ================= MAIN ======================
int main(){

    mt19937 rng(random_device{}());

    L     = sqrt(N / rho);
    halfL = 0.5 * L;

    vector<Particle> P = initialize(N, rng);

    ofstream fout("simulation_data.txt");
    fout << "# L=" << L << "\n";

    build_neighbor_list(P);

    int frame = 0;

    for(int step = 0; step < STEPS; step++){

        if(need_rebuild(P))
            build_neighbor_list(P);

        update(P, rng);

        if(step % SAVE_FREQ == 0)
            save_snapshot(fout, P, frame++);
    }

    fout.close();
    cout << "Simulation complete.\n";
}


