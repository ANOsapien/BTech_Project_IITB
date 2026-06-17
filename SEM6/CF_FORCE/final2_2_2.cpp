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
const int    N =2;
const double v0 =0.2;
const double r_int = 1.0;

const double beta_CF  = 0.0;
const double alpha_CF = 1.0;

const double Dr = 0.0;
const double rho = 0.01;


const double dt = 0.001;
const int    STEPS = 80000;
const int    SAVE_FREQ =10;

// ================= NEIGHBOR LIST ==============
const double r_skin = 0.3;
const double r_cut  = r_int + r_skin;

vector<vector<int>> neighbors;
vector<double> lastx,lasty;

// ================= GLOBALS ====================
double L;

// ================= PARTICLE ==================
struct Particle {
    double x,y;
    double theta_n;
    double vx,vy;
    double fx,fy;
};

// ================= PBC =======================
inline void wrap_position(double &x,double &y){
    if(x<0) x+=L;
    if(x>=L) x-=L;
    if(y<0) y+=L;
    if(y>=L) y-=L;
}

inline void wrap_distance(double &dx,double &dy){
    if(dx> L/2) dx-=L;
    if(dx<-L/2) dx+=L;
    if(dy> L/2) dy-=L;
    if(dy<-L/2) dy+=L;
}

// ================= SNAPSHOT ==================
void save_snapshot(ofstream &fout,const vector<Particle>& P,int frame){

    fout<<"FRAME "<<frame<<"\n";
    for(auto &p:P)
        fout<<p.x<<" "<<p.y<<" "<<p.vx<<" "<<p.vy<<" "
            <<p.theta_n<<" "<<p.fx<<" "<<p.fy<<"\n";
}

// ================= INITIALIZE =================
vector<Particle> initialize(int N, mt19937 &rng){

    vector<Particle> P(N);

    double cy=0.5*L;
    double alpha=0.4;
    double spacing=alpha*2*r_int;

    double x0=0.5*(L-spacing*(N-1));

    for(int i=0;i<N;i++){

        P[i].x=x0+i*spacing;
        P[i].y=cy;
        wrap_position(P[i].x,P[i].y);

        if(i==0) P[i].theta_n=1.57;
        else if(i==N-1) P[i].theta_n=0;
        else P[i].theta_n=0;

        P[i].vx=v0*cos(P[i].theta_n);
        P[i].vy=v0*sin(P[i].theta_n);
    }
    return P;
}


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


// vector<Particle> initialize(int N, mt19937 &rng) {

//     vector<Particle> P(N);

//     double cx = 0.5 * L;
//     double cy = 0.5 * L;

//     // perfect touching ring radius
//     double R = r_int / (2.0 * sin(M_PI / N));

//     uniform_real_distribution<double> noise(-0.01*r_int, 0.01*r_int);

//     for(int i=0;i<N;i++){

//         double phi = 2.0*M_PI*i/N ;

//         // perfect circle
//         double x0 = cx + R*cos(phi);
//         double y0 = cy + R*sin(phi);

//         // small distortion (keeps near touching)
//         P[i].x = x0 + noise(rng);
//         P[i].y = y0 + noise(rng);

//         wrap_position(P[i].x,P[i].y);

//         // tangential polarity
//         P[i].theta_n = phi + M_PI/2;

//         P[i].vx = v0*cos(P[i].theta_n);
//         P[i].vy = v0*sin(P[i].theta_n);
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
//             P[i].theta_n =  1.5;     // +90°
//         else if (i == N-1)
//             P[i].theta_n = 1.5;     // -90°
//         else
//             P[i].theta_n = 0.0;         // straight

//         P[i].vx = v0 * cos(P[i].theta_n);
//         P[i].vy = v0 * sin(P[i].theta_n);
//     }

//     return P;
// }

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


// ================= NEIGHBOR LIST =================
void build_neighbor_list(const vector<Particle>& P){

    int N=P.size();
    neighbors.assign(N,{});

    for(int i=0;i<N;i++)
        for(int j=i+1;j<N;j++){

            double dx=P[j].x-P[i].x;
            double dy=P[j].y-P[i].y;
            wrap_distance(dx,dy);

            if(dx*dx+dy*dy<r_cut*r_cut){
                neighbors[i].push_back(j);
                neighbors[j].push_back(i);
            }
        }

    lastx.resize(N);
    lasty.resize(N);
    for(int i=0;i<N;i++){
        lastx[i]=P[i].x;
        lasty[i]=P[i].y;
    }
}

bool need_rebuild(const vector<Particle>& P){

    for(int i=0;i<P.size();i++){
        double dx=P[i].x-lastx[i];
        double dy=P[i].y-lasty[i];
        wrap_distance(dx,dy);

        if(dx*dx+dy*dy>0.25*r_skin*r_skin)
            return true;
    }
    return false;
}

// ================= FORCE COMPUTE =================
// void compute_forces(const vector<Particle>& P,
//                     vector<pair<double,double>>& forces,
//                     vector<pair<double,double>>& JCF){

//     int N=P.size();

//     for(int i=0;i<N;i++){
//         forces[i]={0,0};
//         JCF[i]={0,0};
//     }

// for(int i=0;i<N;i++)
//     for(int j:neighbors[i]){

//         if(i==j) continue;

//         double dx = P[j].x - P[i].x;
//         double dy = P[j].y - P[i].y;
//         wrap_distance(dx,dy);

//         double d2 = dx*dx + dy*dy;
//         if(d2 < 1e-12 || d2 >= r_int*r_int) continue;

//         double d = sqrt(d2);
//         double dxn = dx/d;
//         double dyn = dy/d;

//         double Fmag = beta_CF * (r_int/d - 1.0);

//         forces[i].first  -= Fmag*dxn;
//         forces[i].second -= Fmag*dyn;
//         forces[j].first  += Fmag*dxn;
//         forces[j].second += Fmag*dyn;

//         // --- NON-RECIPROCAL TORQUE ON i ONLY ---
//         double qjx = cos(P[j].theta_n);
//         double qjy = sin(P[j].theta_n);

//         double w = 0.5*(1 + dxn*qjx + dyn*qjy);

//         JCF[i].first  += w*dxn;
//         JCF[i].second += w*dyn;
//     }
// }

void compute_forces(const vector<Particle>& P,
                    vector<pair<double,double>>& forces,
                    vector<pair<double,double>>& JCF){

    int N=P.size();

    for(int i=0;i<N;i++){
        forces[i]={0,0};
        JCF[i]={0,0};
    }

    for(int i=0;i<N;i++)
        for(int j:neighbors[i]){

            if(j <= i) continue;

            double dx = P[j].x - P[i].x;
            double dy = P[j].y - P[i].y;
            wrap_distance(dx,dy);

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

            // --- NON-RECIPROCAL TORQUE ON i (affected by j's polarity) ---
            double qjx = cos(P[j].theta_n);
            double qjy = sin(P[j].theta_n);

            double wi = 0.5*(1 + dxn*qjx + dyn*qjy);

            JCF[i].first  += wi*dxn;
            JCF[i].second += wi*dyn;

            // --- NON-RECIPROCAL TORQUE ON j (affected by i's polarity) ---
            double qix = cos(P[i].theta_n);
            double qiy = sin(P[i].theta_n);

            double wj = 0.5*(1 - dxn*qix - dyn*qiy);

            JCF[j].first  += wj*(-dxn);
            JCF[j].second += wj*(-dyn);
        }
}





// ================= HEUN UPDATE =================
void update(vector<Particle>& P, mt19937 &rng, int step){

    int N=P.size();

    vector<pair<double,double>> F0(N),J0(N);
    vector<pair<double,double>> F1(N),J1(N);

    compute_forces(P,F0,J0);

    vector<Particle> Ppred=P;

    // predictor
    for(int i=0;i<N;i++){

        double vx0=v0*cos(P[i].theta_n)+F0[i].first;
        double vy0=v0*sin(P[i].theta_n)+F0[i].second;

        Ppred[i].x+=vx0*dt;
        Ppred[i].y+=vy0*dt;
        wrap_position(Ppred[i].x,Ppred[i].y);

        double nx=cos(P[i].theta_n);
        double ny=sin(P[i].theta_n);

        Ppred[i].theta_n+=dt*alpha_CF*(J0[i].first*(-ny)+J0[i].second*nx);
    }

    compute_forces(Ppred,F1,J1);

    // corrector
    for(int i=0;i<N;i++){

        double vx0=v0*cos(P[i].theta_n)+F0[i].first;
        double vy0=v0*sin(P[i].theta_n)+F0[i].second;

        double vx1=v0*cos(Ppred[i].theta_n)+F1[i].first;
        double vy1=v0*sin(Ppred[i].theta_n)+F1[i].second;

        P[i].x+=0.5*(vx0+vx1)*dt;
        P[i].y+=0.5*(vy0+vy1)*dt;
        wrap_position(P[i].x,P[i].y);

        double n0x=cos(P[i].theta_n);
        double n0y=sin(P[i].theta_n);
        double n1x=cos(Ppred[i].theta_n);
        double n1y=sin(Ppred[i].theta_n);

        double d0=alpha_CF*(J0[i].first*(-n0y)+J0[i].second*n0x);
        double d1=alpha_CF*(J1[i].first*(-n1y)+J1[i].second*n1x);

        P[i].theta_n+=0.5*(d0+d1)*dt;

        P[i].vx=vx1;
        P[i].vy=vy1;
        P[i].fx=F1[i].first;
        P[i].fy=F1[i].second;
    }
}

// ================= MAIN ======================
int main(){

    mt19937 rng(random_device{}());

    L=sqrt(N/rho);

    vector<Particle> P=initialize(N,rng);

    ofstream fout("simulation_data.txt");
    fout<<"# L="<<L<<"\n";

    forceout.open("particle1_full_trace.txt");
    forceout<<"# step x y vx vy theta Fx Fy\n";

    build_neighbor_list(P);

    int frame=0;

    for(int step=0;step<STEPS;step++){

        if(need_rebuild(P))
            build_neighbor_list(P);

        update(P,rng,step);

        if(step%SAVE_FREQ==0)
            save_snapshot(fout,P,frame++);
    }

    fout.close();
    forceout.close();

    cout<<"Simulation complete.\n";
}