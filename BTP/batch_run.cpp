#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <string>
using namespace std;

// ================= PARAMETERS =================
const int    N    = 1000;
const double v0   = 1.0;
const double r_int = 1.0;

const double beta_CF  = 1.0;
const double alpha_CF = 1.0;

const double Dr  = 0.1;
const double rho = 1.0;

const double dt        = 0.001;
const int    STEPS     = 500000;
const int    SAVE_FREQ = 10;

// ================= NEIGHBOR LIST ==============
const double r_skin = 0.3;
const double r_cut  = r_int + r_skin;
const double r_cut2 = r_cut * r_cut;
const double r_int2 = r_int * r_int;
const double r_skin2_quarter = 0.25 * r_skin * r_skin;

vector<vector<int>> neighbors;
vector<double> lastx, lasty;

// ================= GLOBALS ====================
double L;
double halfL;

// ================= PARTICLE ==================
struct Particle {
    double x, y;
    double theta_n;
    double cos_theta, sin_theta;
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

inline double wrap_angle(double theta){
    theta = fmod(theta, 2.0*M_PI);
    if(theta < 0) theta += 2.0*M_PI;
    return theta;
}

// ================= RUN FOLDER =================
string make_run_dir(const string &init_name){
    ostringstream oss;
    oss << "run"
        << "_" << init_name
        << "_N"   << N
        << "_rho" << rho
        << "_b"   << beta_CF
        << "_a"   << alpha_CF
        << "_Dr"  << Dr;
    string dir = oss.str();

    string candidate = dir;
    int suffix = 2;
    while(true){
        ifstream test(candidate + "/params.txt");
        if(!test.good()) break;
        ostringstream ss;
        ss << dir << "_" << suffix++;
        candidate = ss.str();
    }

    string cmd = string("mkdir \"") + candidate + string("\"");
    (void)system(cmd.c_str());
    return candidate;
}

void save_params(const string &dir, const string &init_name){
    ofstream f(dir + "/params.txt");
    f << "init_cond = " << init_name  << "\n";
    f << "N         = " << N          << "\n";
    f << "v0        = " << v0         << "\n";
    f << "r_int     = " << r_int      << "\n";
    f << "beta_CF   = " << beta_CF    << "\n";
    f << "alpha_CF  = " << alpha_CF   << "\n";
    f << "Dr        = " << Dr         << "\n";
    f << "rho       = " << rho        << "\n";
    f << "dt        = " << dt         << "\n";
    f << "STEPS     = " << STEPS      << "\n";
    f << "SAVE_FREQ = " << SAVE_FREQ  << "\n";
    f << "L         = " << sqrt(N / rho) << "\n";
    f << "r_skin    = " << r_skin     << "\n";
    f.close();
}

void save_notes_template(const string &dir, const string &init_name){
    ofstream f(dir + "/notes.txt");
    f << "# Notes for run: " << dir << "\n";
    f << "# Initial condition: " << init_name << "\n";
    f << "# Lines starting with # are comments.\n\n";
    f << "Purpose:\n\n";
    f << "Initial condition:\n\n";
    f << "Observations:\n\n";
    f << "Issues / anomalies:\n\n";
    f << "Next steps:\n\n";
    f.close();
}

// ================= SNAPSHOT ==================
void save_snapshot(ofstream &fout, const vector<Particle>& P, int frame){
    fout << "FRAME " << frame << "\n";
    for(auto &p : P)
        fout << p.x  << " " << p.y  << " "
             << p.vx << " " << p.vy << " "
             << p.theta_n << " " << p.fx << " " << p.fy << "\n";
}

// =====================================================
// ================= INIT FUNCTIONS ================
// =====================================================

// 1. Concentric rings — CCW tangential velocity
vector<Particle> init_concentric_rings(mt19937 &rng){
    vector<Particle> P; P.reserve(N);
    double cx = 0.5*L, cy = 0.5*L;
    double compression = 0.5;
    int count = 0, k = 1;
    while(count < N){
        double R = k * r_int * compression;
        int Nk = max(6, (int)round(2.0*M_PI*k));
        for(int i = 0; i < Nk && count < N; i++){
            double phi = 2.0*M_PI*i/Nk;
            Particle p;
            p.x = cx + R*cos(phi);
            p.y = cy + R*sin(phi);
            p.theta_n   = phi + M_PI/2.0;
            p.cos_theta = cos(p.theta_n);
            p.sin_theta = sin(p.theta_n);
            p.vx = v0*p.cos_theta; p.vy = v0*p.sin_theta;
            p.fx = p.fy = 0.0;
            P.push_back(p); count++;
        }
        k++;
    }
    return P;
}

// 2. Lattice — random angles
vector<Particle> init_lattice(mt19937 &rng){
    vector<Particle> P(N);
    int nSide = (int)ceil(sqrt(N));
    double a = L / nSide;
    uniform_real_distribution<double> ang(0.0, 2*M_PI);
    int idx = 0;
    for(int i = 0; i < nSide && idx < N; i++)
        for(int j = 0; j < nSide && idx < N; j++){
            P[idx].x = (i+0.5)*a; P[idx].y = (j+0.5)*a;
            wrap_position(P[idx].x, P[idx].y);
            P[idx].theta_n   = ang(rng);
            P[idx].cos_theta = cos(P[idx].theta_n);
            P[idx].sin_theta = sin(P[idx].theta_n);
            P[idx].vx = v0*P[idx].cos_theta;
            P[idx].vy = v0*P[idx].sin_theta;
            P[idx].fx = P[idx].fy = 0.0;
            idx++;
        }
    return P;
}

// 3. Two clusters — head-on collision
vector<Particle> init_two_clusters(mt19937 &rng){
    vector<Particle> P; P.reserve(N);
    int N_half = N/2;
    double r_cluster = 5.0*r_int;
    double cx_left = 0.2*L, cx_right = 0.8*L, cy = 0.5*L;
    uniform_real_distribution<double> ang(0.0, 2*M_PI);
    uniform_real_distribution<double> rad(0.0, r_cluster);
    for(int i = 0; i < N_half; i++){
        double a = ang(rng), r = sqrt(rad(rng)*rad(rng));
        Particle p;
        p.x = cx_left + r*cos(a); p.y = cy + r*sin(a);
        wrap_position(p.x, p.y);
        p.theta_n = 0.0;
        p.cos_theta = 1.0; p.sin_theta = 0.0;
        p.vx = v0; p.vy = 0.0; p.fx = p.fy = 0.0;
        P.push_back(p);
    }
    for(int i = N_half; i < N; i++){
        double a = ang(rng), r = sqrt(rad(rng)*rad(rng));
        Particle p;
        p.x = cx_right + r*cos(a); p.y = cy + r*sin(a);
        wrap_position(p.x, p.y);
        p.theta_n = M_PI;
        p.cos_theta = -1.0; p.sin_theta = 0.0;
        p.vx = -v0; p.vy = 0.0; p.fx = p.fy = 0.0;
        P.push_back(p);
    }
    return P;
}

// 4. Vertical lines — head-on
vector<Particle> init_vertical_lines(mt19937 &rng){
    vector<Particle> P; P.reserve(N);
    int N_half = N/2;
    double x_left = 0.2*L, x_right = 0.8*L;
    double y_spacing = L / (N_half+1);
    for(int i = 0; i < N_half; i++){
        Particle p;
        p.x = x_left; p.y = (i+1)*y_spacing;
        wrap_position(p.x, p.y);
        p.theta_n = 0.0; p.cos_theta = 1.0; p.sin_theta = 0.0;
        p.vx = v0; p.vy = 0.0; p.fx = p.fy = 0.0;
        P.push_back(p);
    }
    for(int i = N_half; i < N; i++){
        Particle p;
        p.x = x_right; p.y = (i-N_half+1)*y_spacing;
        wrap_position(p.x, p.y);
        p.theta_n = M_PI; p.cos_theta = -1.0; p.sin_theta = 0.0;
        p.vx = -v0; p.vy = 0.0; p.fx = p.fy = 0.0;
        P.push_back(p);
    }
    return P;
}

// 5. Horizontal rows
vector<Particle> init_horizontal_rows(mt19937 &rng){
    vector<Particle> P; P.reserve(N);
    int per_side = N/2;
    int slabs = 6, rows_per_slab = 6;
    double dx = 2.1*r_int, dy = 2.1*r_int;
    double slab_gap = 4*r_int;
    int count = 0;
    for(int s = 0; s < slabs && count < per_side; s++){
        double y0 = (s+1)*slab_gap + s*rows_per_slab*dy;
        for(int r = 0; r < rows_per_slab && count < per_side; r++){
            double y = y0 + r*dy;
            for(int i = 0; count < per_side; i++){
                double x = 0.05*L + i*dx;
                if(x > 0.4*L) break;
                Particle p;
                p.x = x; p.y = y;
                wrap_position(p.x, p.y);
                p.theta_n = 0.0; p.cos_theta = 1.0; p.sin_theta = 0.0;
                p.vx = v0; p.vy = 0.0; p.fx = p.fy = 0.0;
                P.push_back(p); count++;
            }
        }
    }
    for(int s = 0; s < slabs && count < N; s++){
        double y0 = (s+1)*slab_gap + s*rows_per_slab*dy;
        for(int r = 0; r < rows_per_slab && count < N; r++){
            double y = y0 + r*dy;
            for(int i = 0; count < N; i++){
                double x = 0.95*L - i*dx;
                if(x < 0.6*L) break;
                Particle p;
                p.x = x; p.y = y;
                wrap_position(p.x, p.y);
                p.theta_n = M_PI; p.cos_theta = -1.0; p.sin_theta = 0.0;
                p.vx = -v0; p.vy = 0.0; p.fx = p.fy = 0.0;
                P.push_back(p); count++;
            }
        }
    }
    return P;
}

// 6. Two lattice blocks — left vs right
vector<Particle> init_two_rows(mt19937 &rng){
    vector<Particle> P; P.reserve(N);
    int N_half = N/2;
    int nSide = (int)ceil(sqrt(N_half));
    double gap = 0.2*L;
    double width_block = (L - gap)/2.0;
    double a = width_block / nSide;
    int idx = 0;
    for(int i = 0; i < nSide && idx < N_half; i++)
        for(int j = 0; j < nSide && idx < N_half; j++){
            Particle p;
            p.x = (i+0.5)*a; p.y = (j+0.5)*a;
            wrap_position(p.x, p.y);
            p.theta_n = 0.0; p.cos_theta = 1.0; p.sin_theta = 0.0;
            p.vx = v0; p.vy = 0.0; p.fx = p.fy = 0.0;
            P.push_back(p); idx++;
        }
    double x_right = width_block + gap;
    for(int i = 0; i < nSide && idx < N; i++)
        for(int j = 0; j < nSide && idx < N; j++){
            Particle p;
            p.x = x_right + (i+0.5)*a; p.y = (j+0.5)*a;
            wrap_position(p.x, p.y);
            p.theta_n = M_PI; p.cos_theta = -1.0; p.sin_theta = 0.0;
            p.vx = -v0; p.vy = 0.0; p.fx = p.fy = 0.0;
            P.push_back(p); idx++;
        }
    return P;
}

// 7. Four lattice blocks
vector<Particle> init_four_blocks(mt19937 &rng){
    vector<Particle> P; P.reserve(N);
    int per_block = N/4;
    int nSide = (int)ceil(sqrt(per_block));
    double gap_x = 0.2*L, gap_y = 0.2*L;
    double block_w = (L - gap_x)/2.0, block_h = (L - gap_y)/2.0;
    double a = min(block_w, block_h) / nSide;
    double xs[4] = {0.0, block_w+gap_x, 0.0, block_w+gap_x};
    double ys[4] = {0.0, 0.0, block_h+gap_y, block_h+gap_y};
    double thetas[4] = {0.0, M_PI, 0.0, M_PI};
    double vxs[4] = {v0, -v0, v0, -v0};
    for(int b = 0; b < 4; b++){
        int count = 0;
        for(int i = 0; i < nSide && count < per_block; i++)
            for(int j = 0; j < nSide && count < per_block; j++){
                Particle p;
                p.x = xs[b]+(i+0.5)*a; p.y = ys[b]+(j+0.5)*a;
                wrap_position(p.x, p.y);
                p.theta_n = thetas[b];
                p.cos_theta = cos(p.theta_n); p.sin_theta = sin(p.theta_n);
                p.vx = vxs[b]; p.vy = 0.0; p.fx = p.fy = 0.0;
                P.push_back(p); count++;
            }
    }
    return P;
}

// 8. Hexagonal lattice
vector<Particle> init_hex_lattice(mt19937 &rng){
    vector<Particle> P; P.reserve(N);
    int sites = N/6;
    double a = L/sqrt(sites);
    double dy = a*sqrt(3.0)/2.0;
    double r_hex = 0.4*r_int;
    uniform_real_distribution<double> ang(0.0, 2*M_PI);
    int idx = 0, row = 0;
    for(double y = dy/2; y < L && idx < N; y += dy, row++){
        double x_offset = (row%2) ? a/2 : 0;
        for(double x = a/2+x_offset; x < L && idx < N; x += a){
            double theta_site = ang(rng);
            for(int k = 0; k < 6 && idx < N; k++){
                double phi = 2*M_PI*k/6.0;
                Particle p;
                p.x = x + r_hex*cos(phi); p.y = y + r_hex*sin(phi);
                wrap_position(p.x, p.y);
                p.theta_n   = theta_site;
                p.cos_theta = cos(p.theta_n); p.sin_theta = sin(p.theta_n);
                p.vx = v0*p.cos_theta; p.vy = v0*p.sin_theta;
                p.fx = p.fy = 0.0;
                P.push_back(p); idx++;
            }
        }
    }
    return P;
}

// 9. Shear flow
vector<Particle> init_shear_flow(mt19937 &rng){
    vector<Particle> P(N);
    int nSide = (int)ceil(sqrt(N));
    double a = L/nSide;
    double gamma = 1.0;
    int idx = 0;
    for(int i = 0; i < nSide && idx < N; i++)
        for(int j = 0; j < nSide && idx < N; j++){
            P[idx].x = (i+0.5)*a; P[idx].y = (j+0.5)*a;
            wrap_position(P[idx].x, P[idx].y);
            double vx = gamma*(P[idx].y - 0.5*L), vy = 0.0;
            P[idx].vx = vx; P[idx].vy = vy;
            P[idx].theta_n   = atan2(vy, vx);
            P[idx].cos_theta = cos(P[idx].theta_n);
            P[idx].sin_theta = sin(P[idx].theta_n);
            P[idx].fx = P[idx].fy = 0.0;
            idx++;
        }
    return P;
}

// 10. Four lines
vector<Particle> init_four_lines(mt19937 &rng){
    vector<Particle> P; P.reserve(N);
    double spacing_x = 0.5*r_int, spacing_y = 1.5*r_int;
    double offset_x = 0.25*L, offset_y = 0.15*L;
    int per_block = N/4;
    int nSide = (int)ceil(sqrt((double)per_block));
    double cx = 0.5*L, cy = 0.5*L;
    double bx[4]    = {cx-offset_x, cx+offset_x, cx-offset_x, cx+offset_x};
    double by[4]    = {cy+offset_y, cy+offset_y, cy-offset_y, cy-offset_y};
    double bvx[4]   = {v0, -v0, v0, -v0};
    double btheta[4]= {0.0, M_PI, 0.0, M_PI};
    for(int b = 0; b < 4; b++){
        int count = 0;
        double x0 = bx[b] - 0.5*(nSide-1)*spacing_x;
        double y0 = by[b] - 0.5*(nSide-1)*spacing_y;
        for(int i = 0; i < nSide && count < per_block; i++)
            for(int j = 0; j < nSide && count < per_block; j++){
                Particle p;
                p.x = x0+i*spacing_x; p.y = y0+j*spacing_y;
                wrap_position(p.x, p.y);
                p.theta_n   = btheta[b];
                p.cos_theta = cos(p.theta_n); p.sin_theta = sin(p.theta_n);
                p.vx = bvx[b]; p.vy = 0.0; p.fx = p.fy = 0.0;
                P.push_back(p); count++;
            }
    }
    return P;
}

// 11. Out velocity — outermost ring moving tangentially
vector<Particle> init_out_velocity(mt19937 &rng){
    vector<Particle> P; P.reserve(N);
    double cx = 0.5*L, cy = 0.5*L;
    double compression = 0.5;
    double omega = 1.0;
    int count = 0, k = 1;
    int last_ring_start = 0;
    while(count < N){
        double R = k*r_int*compression;
        int Nk = max(6, (int)round(2.0*M_PI*k));
        last_ring_start = count;
        for(int i = 0; i < Nk && count < N; i++){
            double phi = 2.0*M_PI*i/Nk;
            Particle p;
            p.x = cx + R*cos(phi); p.y = cy + R*sin(phi);
            p.theta_n = phi + M_PI/2.0;
            p.cos_theta = cos(p.theta_n); p.sin_theta = sin(p.theta_n);
            p.vx = 0.0; p.vy = 0.0; p.fx = p.fy = 0.0;
            P.push_back(p); count++;
        }
        k++;
    }
    // give outermost ring tangential velocity
    for(int j = last_ring_start; j < count; j++){
        double dx = P[j].x - cx, dy = P[j].y - cy;
        double phi = atan2(dy, dx);
        double R_ring = sqrt(dx*dx + dy*dy);
        double v = omega*R_ring;
        P[j].vx = -v*sin(phi);
        P[j].vy =  v*cos(phi);
    }
    return P;
}

// 12. In velocity — innermost ring moving tangentially
vector<Particle> init_in_velocity(mt19937 &rng){
    vector<Particle> P; P.reserve(N);
    double cx = 0.5*L, cy = 0.5*L;
    double compression = 1.1;
    double omega = 1.0;
    int count = 0, k = 1;
    while(count < N){
        double R = k*r_int*compression;
        int Nk = max(6, (int)round(2.0*M_PI*k));
        for(int i = 0; i < Nk && count < N; i++){
            double phi = 2.0*M_PI*i/Nk;
            Particle p;
            p.x = cx + R*cos(phi); p.y = cy + R*sin(phi);
            p.theta_n = phi + M_PI/2.0;
            p.cos_theta = cos(p.theta_n); p.sin_theta = sin(p.theta_n);
            if(k == 1){
                double v = omega*R;
                p.vx = -v*sin(phi); p.vy = v*cos(phi);
            } else {
                p.vx = 0.0; p.vy = 0.0;
            }
            p.fx = p.fy = 0.0;
            P.push_back(p); count++;
        }
        k++;
    }
    return P;
}

// 13. Single flock — all aligned in one direction
vector<Particle> init_single_flock(mt19937 &rng){
    vector<Particle> P(N);
    int nSide = (int)ceil(sqrt(N));
    double a = L/nSide;
    int idx = 0;
    for(int i = 0; i < nSide && idx < N; i++)
        for(int j = 0; j < nSide && idx < N; j++){
            P[idx].x = (i+0.5)*a; P[idx].y = (j+0.5)*a;
            wrap_position(P[idx].x, P[idx].y);
            P[idx].theta_n   = 0.0;   // all pointing right
            P[idx].cos_theta = 1.0; P[idx].sin_theta = 0.0;
            P[idx].vx = v0; P[idx].vy = 0.0; P[idx].fx = P[idx].fy = 0.0;
            idx++;
        }
    return P;
}

// 14. Vortex — spinning disk, angular velocity proportional to radius
vector<Particle> init_vortex(mt19937 &rng){
    vector<Particle> P; P.reserve(N);
    double cx = 0.5*L, cy = 0.5*L;
    double compression = 0.5;
    double omega = 1.0;
    int count = 0, k = 1;
    while(count < N){
        double R = k*r_int*compression;
        int Nk = max(6, (int)round(2.0*M_PI*k));
        for(int i = 0; i < Nk && count < N; i++){
            double phi = 2.0*M_PI*i/Nk;
            Particle p;
            p.x = cx + R*cos(phi); p.y = cy + R*sin(phi);
            p.theta_n   = phi + M_PI/2.0;   // tangential
            p.cos_theta = cos(p.theta_n); p.sin_theta = sin(p.theta_n);
            double v = omega*R;              // speed proportional to radius
            p.vx = -v*sin(phi); p.vy = v*cos(phi);
            p.fx = p.fy = 0.0;
            P.push_back(p); count++;
        }
        k++;
    }
    return P;
}

// 15. Anti-vortex — two halves spinning in opposite directions
vector<Particle> init_anti_vortex(mt19937 &rng){
    vector<Particle> P; P.reserve(N);
    double cx = 0.5*L, cy = 0.5*L;
    double compression = 0.5;
    double omega = 1.0;
    int count = 0, k = 1;
    while(count < N){
        double R = k*r_int*compression;
        int Nk = max(6, (int)round(2.0*M_PI*k));
        for(int i = 0; i < Nk && count < N; i++){
            double phi = 2.0*M_PI*i/Nk;
            Particle p;
            p.x = cx + R*cos(phi); p.y = cy + R*sin(phi);
            // top half CCW, bottom half CW
            double sign = (p.y > cy) ? 1.0 : -1.0;
            p.theta_n   = phi + sign*M_PI/2.0;
            p.cos_theta = cos(p.theta_n); p.sin_theta = sin(p.theta_n);
            double v = omega*R;
            p.vx = -sign*v*sin(phi); p.vy = sign*v*cos(phi);
            p.fx = p.fy = 0.0;
            P.push_back(p); count++;
        }
        k++;
    }
    return P;
}

// ================= DISPATCHER =================
vector<Particle> initialize(const string &name, mt19937 &rng){
    if(name == "rings")         return init_concentric_rings(rng);
    if(name == "lattice")       return init_lattice(rng);
    if(name == "two_clusters")  return init_two_clusters(rng);
    if(name == "vert_lines")    return init_vertical_lines(rng);
    if(name == "horiz_rows")    return init_horizontal_rows(rng);
    if(name == "two_rows")      return init_two_rows(rng);
    if(name == "four_blocks")   return init_four_blocks(rng);
    if(name == "hex_lattice")   return init_hex_lattice(rng);
    if(name == "shear_flow")    return init_shear_flow(rng);
    if(name == "four_lines")    return init_four_lines(rng);
    if(name == "out_velocity")  return init_out_velocity(rng);
    if(name == "in_velocity")   return init_in_velocity(rng);
    if(name == "single_flock")  return init_single_flock(rng);
    if(name == "vortex")        return init_vortex(rng);
    if(name == "anti_vortex")   return init_anti_vortex(rng);
    cerr << "Unknown init condition: " << name << "\n";
    exit(1);
}

// ================= NEIGHBOR LIST =================
void build_neighbor_list(const vector<Particle>& P){
    int n = P.size();
    neighbors.assign(n, {});
    for(int i = 0; i < n; i++)
        for(int j = i+1; j < n; j++){
            double dx = P[j].x - P[i].x, dy = P[j].y - P[i].y;
            wrap_distance(dx, dy);
            if(dx*dx + dy*dy < r_cut2){
                neighbors[i].push_back(j);
                neighbors[j].push_back(i);
            }
        }
    lastx.resize(n); lasty.resize(n);
    for(int i = 0; i < n; i++){ lastx[i] = P[i].x; lasty[i] = P[i].y; }
}

bool need_rebuild(const vector<Particle>& P){
    for(int i = 0; i < (int)P.size(); i++){
        double dx = P[i].x - lastx[i], dy = P[i].y - lasty[i];
        wrap_distance(dx, dy);
        if(dx*dx + dy*dy > r_skin2_quarter) return true;
    }
    return false;
}

// ================= FORCE COMPUTE =================
void compute_forces(const vector<Particle>& P,
                    vector<pair<double,double>>& forces,
                    vector<pair<double,double>>& JCF){
    int n = P.size();
    for(int i = 0; i < n; i++){ forces[i] = {0,0}; JCF[i] = {0,0}; }
    for(int i = 0; i < n; i++){
        for(int j : neighbors[i]){
            if(j <= i) continue;
            double dx = P[j].x-P[i].x, dy = P[j].y-P[i].y;
            wrap_distance(dx, dy);
            double d2 = dx*dx + dy*dy;
            if(d2 < 1e-12 || d2 >= r_int2) continue;
            double inv_d = 1.0/sqrt(d2);
            double dxn = dx*inv_d, dyn = dy*inv_d;
            if(beta_CF != 0.0){
                double Fmag = beta_CF*(r_int*inv_d - 1.0);
                forces[i].first  -= Fmag*dxn; forces[i].second -= Fmag*dyn;
                forces[j].first  += Fmag*dxn; forces[j].second += Fmag*dyn;
            }
            if(alpha_CF != 0.0){
                double dot_j = dxn*P[j].cos_theta + dyn*P[j].sin_theta;
                double wi = 0.5*(1.0+dot_j);
                JCF[i].first += wi*dxn; JCF[i].second += wi*dyn;
                double dot_i = dxn*P[i].cos_theta + dyn*P[i].sin_theta;
                double wj = 0.5*(1.0-dot_i);
                JCF[j].first -= wj*dxn; JCF[j].second -= wj*dyn;
            }
        }
    }
}

// ================= HEUN UPDATE =================
void update(vector<Particle>& P, mt19937 &rng){
    int n = P.size();
    static vector<pair<double,double>> F0, J0, F1, J1;
    static vector<Particle> Ppred;
    F0.resize(n); J0.resize(n); F1.resize(n); J1.resize(n); Ppred.resize(n);
    static normal_distribution<double> gauss(0.0, 1.0);
    double noise_amp = sqrt(2.0*Dr*dt);
    compute_forces(P, F0, J0);
    Ppred = P;
    for(int i = 0; i < n; i++){
        double vx0 = v0*P[i].cos_theta + F0[i].first;
        double vy0 = v0*P[i].sin_theta + F0[i].second;
        Ppred[i].x += vx0*dt; Ppred[i].y += vy0*dt;
        wrap_position(Ppred[i].x, Ppred[i].y);
        double dtheta = dt*alpha_CF*(J0[i].first*(-P[i].sin_theta) + J0[i].second*P[i].cos_theta);
        Ppred[i].theta_n += dtheta;
        Ppred[i].cos_theta = cos(Ppred[i].theta_n);
        Ppred[i].sin_theta = sin(Ppred[i].theta_n);
    }
    compute_forces(Ppred, F1, J1);
    for(int i = 0; i < n; i++){
        double vx0 = v0*P[i].cos_theta    + F0[i].first;
        double vy0 = v0*P[i].sin_theta    + F0[i].second;
        double vx1 = v0*Ppred[i].cos_theta + F1[i].first;
        double vy1 = v0*Ppred[i].sin_theta + F1[i].second;
        P[i].x += 0.5*(vx0+vx1)*dt; P[i].y += 0.5*(vy0+vy1)*dt;
        wrap_position(P[i].x, P[i].y);
        double d0 = alpha_CF*(J0[i].first*(-P[i].sin_theta)    + J0[i].second*P[i].cos_theta);
        double d1 = alpha_CF*(J1[i].first*(-Ppred[i].sin_theta) + J1[i].second*Ppred[i].cos_theta);
        P[i].theta_n += 0.5*(d0+d1)*dt + noise_amp*gauss(rng);
        P[i].theta_n  = wrap_angle(P[i].theta_n);
        P[i].cos_theta = cos(P[i].theta_n); P[i].sin_theta = sin(P[i].theta_n);
        P[i].vx = vx1; P[i].vy = vy1;
        P[i].fx = F1[i].first; P[i].fy = F1[i].second;
    }
}

// ================= MAIN ======================
int main(int argc, char* argv[]){

    // default init or pick from command line
    string init_name = "rings";
    if(argc >= 2) init_name = argv[1];

    mt19937 rng(random_device{}());

    L     = sqrt(N / rho);
    halfL = 0.5 * L;

    string run_dir = make_run_dir(init_name);
    save_params(run_dir, init_name);
    save_notes_template(run_dir, init_name);
    cout << "Init: " << init_name << "  |  Run folder: " << run_dir << "\n";

    vector<Particle> P = initialize(init_name, rng);

    ofstream fout(run_dir + "/simulation.txt");
    fout << "# L=" << L << "\n";

    build_neighbor_list(P);

    int frame = 0;
    save_snapshot(fout, P, frame++);

    for(int step = 0; step < STEPS; step++){
        if(need_rebuild(P)) build_neighbor_list(P);
        update(P, rng);
        if(step % SAVE_FREQ == 0){
            save_snapshot(fout, P, frame++);
            if(frame % 1000 == 0)
                cout << "  frame " << frame << " / " << STEPS/SAVE_FREQ << "\r" << flush;
        }
    }

    fout.close();
    cout << "\nSimulation complete -> " << run_dir << "/\n";

    cout << "Converting to .xyz...\n";
    string conv_cmd = "python convert.py \"" + run_dir + "\"";
    (void)system(conv_cmd.c_str());

    return 0;
}