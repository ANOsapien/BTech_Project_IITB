#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <iomanip>
using namespace std;

// ================= PARAMETERS =================
const int    N    = 900;
const double v0   = 1.0;
const double r_int = 1.0;

const double beta_CF  = 1.0;
const double alpha_CF = 1.0;

// Kuramoto coupling strength and range
const double tau_kur  = 1;         // tau in the Kuramoto equation
const double r_kur    = r_int; // interaction range = 2r

const double alpha_new = 1.0;

const double Dr  = 0;
const double rho = 0.01;

const double dt        = 0.001;
const int    STEPS     = 8000;
const int    SAVE_FREQ = 10;

// ================= NEIGHBOR LIST ==============
const double r_skin = 0.3;

// --- existing list (range r_int) ---
const double r_cut  = r_int + r_skin;
const double r_cut2 = r_cut * r_cut;
const double r_int2 = r_int * r_int;

// --- new Kuramoto list (range 2*r_int) ---
const double r_cut_kur  = r_kur + r_skin;
const double r_cut_kur2 = r_cut_kur * r_cut_kur;
const double r_kur2     = r_kur * r_kur;

const double r_skin2_quarter = 0.25 * r_skin * r_skin;

vector<vector<int>> neighbors;      // within r_int
vector<vector<int>> neighbors_kur;  // within 2*r_int
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
string make_run_dir(){
    ostringstream oss;
    oss << "run"
        << "_N"   << N
        << "_rho" << rho
        << "_b"   << beta_CF
        << "_a"   << alpha_CF
        << "_Dr"  << Dr
        << "_dt"  << dt;
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

    string cmd = string("mkdir \"") + candidate + string("\""); (void)system(cmd.c_str());
    return candidate;
}

void save_params(const string &dir){
    ofstream f(dir + "/params.txt");
    f << "N         = " << N         << "\n";
    f << "v0        = " << v0        << "\n";
    f << "r_int     = " << r_int     << "\n";
    f << "beta_CF   = " << beta_CF   << "\n";
    f << "alpha_CF  = " << alpha_CF  << "\n";
    f << "tau_kur   = " << tau_kur   << "\n";
    f << "r_kur     = " << r_kur     << "\n";
    f << "alpha_new = " << alpha_new << "\n";
    f << "Dr        = " << Dr        << "\n";
    f << "rho       = " << rho       << "\n";
    f << "dt        = " << dt        << "\n";
    f << "STEPS     = " << STEPS     << "\n";
    f << "SAVE_FREQ = " << SAVE_FREQ << "\n";
    f << "L         = " << sqrt(N / rho) << "\n";
    f << "r_skin    = " << r_skin    << "\n";
    f.close();
}

void save_notes_template(const string &dir){
    ofstream f(dir + "/notes.txt");
    f << "# Notes for run: " << dir << "\n";
    f << "# Lines starting with # are comments.\n";
    f << "\n";
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

// ================= INITIALIZE =================
// vector<Particle> initialize(int N, mt19937 &rng){

//     vector<Particle> P(N);

//     int nSide = ceil(sqrt(N));
//     double a = L / nSide;

//     double gamma = 1.0;

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


// concentric ring with velocity
vector<Particle> initialize(int N, mt19937 &rng){

    vector<Particle> P;
    P.reserve(N);

    double cx = 0.5 * L;
    double cy = 0.5 * L;

    double compression = 0.5;   // <1 gives overlap

    int count = 0;
    int k = 1;

    while(count < N){

        double R = k * r_int * compression;

        int Nk = int(round(2.0 * M_PI * k));
        if(Nk < 6) Nk = 6;

        for(int i = 0; i < Nk && count < N; i++){

            double phi = 2.0 * M_PI * i / Nk;

            Particle p;

            p.x = cx + R * cos(phi);
            p.y = cy + R * sin(phi);
            wrap_position(p.x, p.y);

            // tangential direction = phi + pi/2 (CCW)
            p.theta_n   = phi + M_PI / 2.0;
            p.cos_theta = cos(p.theta_n);
            p.sin_theta = sin(p.theta_n);

            // velocity aligned with theta_n — tangential to ring
            p.vx = v0 * p.cos_theta;   // = -v0 * sin(phi)
            p.vy = v0 * p.sin_theta;   // =  v0 * cos(phi)

            p.fx = p.fy = 0.0;

            P.push_back(p);
            count++;
        }
        k++;
    }

    return P;
}



// vector<Particle> initialize(int N, mt19937 &rng){

//     vector<Particle> P;
//     P.reserve(N);

//     double compression = 0.95;
//     double gap         = compression * 2.0 * r_int;

//     double total_len = (N - 1) * gap;
//     double x0        = 0.5 * L - 0.5 * total_len;
//     double y0        = 0.5 * L;

//     // --- tune these directly here ---
//     double theta_first = M_PI/2;
//     double theta_last  = 0.0;
//     double v_first     = 1.0;
//     double v_last      = 1.0;

//     for(int i = 0; i < N; i++){

//         Particle p;
//         p.x = x0 + i * gap;
//         p.y = y0;
//         wrap_position(p.x, p.y);

//         double theta = (i == 0)     ? theta_first :
//                        (i == N - 1) ? theta_last  :
//                                       M_PI / 2.0;

//         double speed = (i == 0)     ? v_first :
//                        (i == N - 1) ? v_last  :
//                                       v0;

//         p.theta_n   = theta;
//         p.cos_theta = cos(theta);
//         p.sin_theta = sin(theta);
//         p.vx        = speed * p.cos_theta;
//         p.vy        = speed * p.sin_theta;
//         p.fx = p.fy = 0.0;

//         P.push_back(p);
//     }

//     return P;
// }

// vector<Particle> initialize(int N, mt19937 &rng){

//     vector<Particle> P;
//     P.reserve(N);

//     double compression = 0.5;
//     double gap         = compression * r_int;

//     double total_len = (N - 1) * gap;
//     double x0        = 0.5 * L - 0.5 * total_len;
//     double y0        = 0.5 * L;

//     double theta_first  = M_PI/2;
//     double theta_middle = 0;
//     double theta_last   = M_PI/2;

//     for(int i = 0; i < N; i++){

//         Particle p;
//         p.x = x0 + i * gap;
//         p.y = y0;
//         wrap_position(p.x, p.y);

//         double theta = (i == 0)     ? theta_first  :
//                        (i == N - 1) ? theta_last   :
//                                       theta_middle;

//         p.theta_n   = theta;
//         p.cos_theta = cos(theta);
//         p.sin_theta = sin(theta);
//         p.vx        = v0 * p.cos_theta;
//         p.vy        = v0 * p.sin_theta;
//         p.fx = p.fy = 0.0;

//         P.push_back(p);
//     }

//     return P;
// }


// ================= NEIGHBOR LIST =================
void build_neighbor_list(const vector<Particle>& P){

    int n = P.size();
    neighbors.assign(n, {});
    neighbors_kur.assign(n, {});   // <-- new

    for(int i = 0; i < n; i++)
        for(int j = i + 1; j < n; j++){
            double dx = P[j].x - P[i].x;
            double dy = P[j].y - P[i].y;
            wrap_distance(dx, dy);

            double d2 = dx*dx + dy*dy;

            if(d2 < r_cut2){
                neighbors[i].push_back(j);
                neighbors[j].push_back(i);
            }

            if(d2 < r_cut_kur2){
                neighbors_kur[i].push_back(j);
                neighbors_kur[j].push_back(i);
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
                    vector<pair<double,double>>& JCF,
                    vector<double>& JKur,
                    vector<pair<double,double>>& JNew){

    int n = P.size();

    for(int i = 0; i < n; i++){
        forces[i] = {0.0, 0.0};
        JCF[i]    = {0.0, 0.0};
        JKur[i]   = 0.0;
        JNew[i]   = {0.0, 0.0};
    }

    for(int i = 0; i < n; i++){
        for(int j : neighbors[i]){

            if(j <= i) continue;

            double dx = P[j].x - P[i].x;
            double dy = P[j].y - P[i].y;
            wrap_distance(dx, dy);

            double d2 = dx*dx + dy*dy;
            if(d2 < 1e-12 || d2 >= r_int2) continue;

            double inv_d = 1.0 / sqrt(d2);
            double dxn   = dx * inv_d;
            double dyn   = dy * inv_d;

            if(beta_CF != 0.0){
                double Fmag = beta_CF * (r_int * inv_d - 1.0);
                forces[i].first  -= Fmag * dxn;
                forces[i].second -= Fmag * dyn;
                forces[j].first  += Fmag * dxn;
                forces[j].second += Fmag * dyn;
            }

            if(alpha_CF != 0.0){
                double dot_j = dxn * P[j].cos_theta + dyn * P[j].sin_theta;
                double wi    = 0.5 * (1.0 + dot_j);
                JCF[i].first  += wi * dxn;
                JCF[i].second += wi * dyn;

                double dot_i = dxn * P[i].cos_theta + dyn * P[i].sin_theta;
                double wj    = 0.5 * (1.0 - dot_i);
                JCF[j].first  -= wj * dxn;
                JCF[j].second -= wj * dyn;
            }
        }
    }

    for(int i = 0; i < n; i++){
        for(int j : neighbors_kur[i]){

            if(j <= i) continue;

            double dx = P[j].x - P[i].x;
            double dy = P[j].y - P[i].y;
            wrap_distance(dx, dy);

            double d2 = dx*dx + dy*dy;
            if(d2 < 1e-12 || d2 >= r_kur2) continue;

            double inv_d = 1.0 / sqrt(d2);
            double dxn   = dx * inv_d;
            double dyn   = dy * inv_d;

            double s = sin(P[j].theta_n - P[i].theta_n);
            JKur[i] += s;
            JKur[j] -= s;

            if(alpha_new != 0.0){
                // double dot_i    = dxn * P[i].cos_theta + dyn * P[i].sin_theta;
                // double dot_j    = dxn * P[j].cos_theta + dyn * P[j].sin_theta;
                // double dot_ninj = P[i].cos_theta * P[j].cos_theta
                //                 + P[i].sin_theta * P[j].sin_theta;

                // double shared = 1.0 - dot_ninj;

                // double wi = 0.25 * (1.0 - dot_i) * shared;
                // JNew[i].first  += wi * dxn;
                // JNew[i].second += wi * dyn;

                // double wj = 0.25 * (1.0 + dot_j) * shared;
                // JNew[j].first  -= wj * dxn;
                // JNew[j].second -= wj * dyn;

                double dot_i = dxn * P[i].cos_theta + dyn * P[i].sin_theta;
                double dot_j = dxn * P[j].cos_theta + dyn * P[j].sin_theta;
                double wi    = 0.5 * (1.0 + dot_i);
                JNew[i].first  += wi * dxn;
                JNew[i].second += wi * dyn;

                double wj    = 0.5 * (1.0 - dot_j);
                JNew[j].first  -= wj * dxn;
                JNew[j].second -= wj * dyn;
            }
        }
    }
}

void update(vector<Particle>& P, mt19937 &rng){

    int n = P.size();

    static vector<pair<double,double>> F0, J0, F1, J1;
    static vector<double> JK0, JK1;
    static vector<pair<double,double>> JN0, JN1;
    static vector<Particle> Ppred;

    F0.resize(n); J0.resize(n); JK0.resize(n); JN0.resize(n);
    F1.resize(n); J1.resize(n); JK1.resize(n); JN1.resize(n);
    Ppred.resize(n);

    static normal_distribution<double> gauss(0.0,1.0);
    double noise_amp = sqrt(2.0 * Dr * dt);

    compute_forces(P, F0, J0, JK0, JN0);

    Ppred = P;

    for(int i = 0; i < n; i++){
        double vx0 = v0 * P[i].cos_theta + F0[i].first;
        double vy0 = v0 * P[i].sin_theta + F0[i].second;

        Ppred[i].x += vx0 * dt;
        Ppred[i].y += vy0 * dt;
        wrap_position(Ppred[i].x, Ppred[i].y);

        double dtheta = dt * alpha_CF * (J0[i].first  * (-P[i].sin_theta)
                                       + J0[i].second *   P[i].cos_theta)
                      + dt * (1.0 * tau_kur) * JK0[i]
                      + dt * alpha_new       * (JN0[i].first  * (-P[i].sin_theta)
                                              + JN0[i].second *   P[i].cos_theta);

        Ppred[i].theta_n += dtheta;
        Ppred[i].cos_theta = cos(Ppred[i].theta_n);
        Ppred[i].sin_theta = sin(Ppred[i].theta_n);
    }

    compute_forces(Ppred, F1, J1, JK1, JN1);

    for(int i = 0; i < n; i++){
        double vx0 = v0 * P[i].cos_theta     + F0[i].first;
        double vy0 = v0 * P[i].sin_theta     + F0[i].second;
        double vx1 = v0 * Ppred[i].cos_theta + F1[i].first;
        double vy1 = v0 * Ppred[i].sin_theta + F1[i].second;

        P[i].x += 0.5 * (vx0 + vx1) * dt;
        P[i].y += 0.5 * (vy0 + vy1) * dt;
        wrap_position(P[i].x, P[i].y);

        double d0 = alpha_CF * (J0[i].first  * (-P[i].sin_theta)
                              + J0[i].second *   P[i].cos_theta)
                  + (1.0 * tau_kur) * JK0[i]
                  + alpha_new       * (JN0[i].first  * (-P[i].sin_theta)
                                     + JN0[i].second *   P[i].cos_theta);

        double d1 = alpha_CF * (J1[i].first  * (-Ppred[i].sin_theta)
                              + J1[i].second *   Ppred[i].cos_theta)
                  + (1.0 * tau_kur) * JK1[i]
                  + alpha_new       * (JN1[i].first  * (-Ppred[i].sin_theta)
                                     + JN1[i].second *   Ppred[i].cos_theta);

        double noise = noise_amp * gauss(rng);
        P[i].theta_n += 0.5 * (d0 + d1) * dt + noise;
        P[i].theta_n  = wrap_angle(P[i].theta_n);

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

    string run_dir = make_run_dir();
    save_params(run_dir);
    save_notes_template(run_dir);
    cout << "Run folder: " << run_dir << "\n";

    vector<Particle> P = initialize(N, rng);

    ofstream fout(run_dir + "/simulation.txt");
    fout << "# L=" << L << "\n";

    build_neighbor_list(P);

    int frame = 0;

    save_snapshot(fout, P, frame++);

    for(int step = 0; step < STEPS; step++){

        if(need_rebuild(P))
            build_neighbor_list(P);

        update(P, rng);

        if(step % SAVE_FREQ == 0)
            save_snapshot(fout, P, frame++);
    }

    fout.close();
    cout << "Simulation complete. Data saved to: " << run_dir << "/\n";

    // Auto-convert simulation.txt -> simulation.xyz
    cout << "Converting to .xyz...\n";
    string conv_cmd = "python convert.py \"" + run_dir + "\"";
    (void)system(conv_cmd.c_str());
}