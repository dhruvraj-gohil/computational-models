#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>

using namespace std;

vector<double> thomas(const vector<double>& a, const vector<double>& b, const vector<double>& c, const vector<double>& d) {
    int n=b.size();
    vector<double> bp=b;
    vector<double> cp=c;
    vector<double> dp=d;

    for (int i=1; i<n; ++i) {
        double m=a[i]/bp[i-1];
        bp[i]-=m*cp[i-1];
        dp[i]-=m*dp[i-1];
    }
    vector<double> x(n);
    x[n-1]=dp[n-1]/bp[n-1];

    for (int i=n-2; i >= 0; --i)
        x[i]=(dp[i]-cp[i]*x[i+1])/bp[i];
    return x;
}

double sim_3_8(const vector<double>& f, double h) {
    int N=f.size()-1;
    if (N%3 != 0) {
        cerr << "Simpson 3/8 requires N%3==0" << endl;
        exit(1);
    }
    double s0=f[0]+f[N];
    double s1=0.0;
    double s2=0.0;

    for (int i=1; i<N; i++) {
        if (i % 3 == 0) s2 += f[i];
        else s1 += f[i];
    }
    return (3.0*h/8.0)*(s0+3.0*s1+2.0*s2);
}

double richardson(const vector<double>& r, const vector<double>& C, double k) {
    const double PI=acos(-1.0);
    int M=r.size();
    int N=M-1;

    if (N%3 != 0) {
        cerr << "Grid must satisfy N%3==0" << endl;
        exit(1);
    }
    double h=r[1]-r[0];

    vector<double> f(M);
    for (int i=0; i<M; i++)
        f[i]=k*C[i]*4.0*PI*r[i]*r[i];

    double I_h=sim_3_8(f, h);

    int M2=2*N+1;
    vector<double> r2(M2), C2(M2), f2(M2);

    for (int j=0; j<M2; j++){
        r2[j]=j*(h/2.0);
    }
    for (int j=0; j<M2; j++) {
        if (j % 2 == 0) C2[j]=C[j/2];
        else C2[j]=0.5*(C[j/2]+C[j/2+1]);
    }

    for (int j=0; j<M2; j++){
        f2[j]=k*C2[j]*4.0*PI*r2[j]*r2[j];
    }
    double I_h2=sim_3_8(f2, h/2.0);
    return I_h2+(I_h2-I_h)/15.0;
}

int main() {
    double R=1e-3;
    double De=1e-6;
    double k=1.0;
    double Cs=1.0;
    double eps=1.0;
    int P=201;
    if ((P-1) % 3 != 0)
        P=((P-1)/3)*3+1;

    cout << "Nodes=" << P << ' ' << endl;

    double dr=R/(P-1);
    vector<double> r(P);
    for (int i=0; i<P; i++) r[i]=i*dr;

    double t_final=5.0;
    double dt=0.01;
    int tp=static_cast<int>(t_final/dt+1e-9); 

    vector<bool> save_steps(tp+1, false);
    double save_points[]={0.0, 0.1, 0.5, 1.0, 2.0, 5.0};
    for (double T : save_points) {
        int idx=static_cast<int>(T/dt+1e-9);
        if (idx >= 0 && idx <= tp) save_steps[idx]=true;
    }

    vector<vector<double>> saved_profiles1;
    vector<double> C(P, 0.0);
    C[P-1]=Cs;
    vector<double> a(P, 0.0), b(P, 0.0), c(P, 0.0);

    //Boundary condition
    b[0]=eps/dt+2*De/(dr*dr)+k;
    c[0]=-2*De/(dr*dr);

    for (int i=1; i<P-1; i++) {
        double ri=r[i];
        a[i]=-De/(dr*dr)-De/(ri*dr);
        b[i]=eps/dt+2*De/(dr*dr)+k;
        c[i]=-De/(dr*dr)+De/(ri*dr);
    }

    b[P-1]=1.0;

    ofstream eta2_out("eta2.csv");
    eta2_out << "time,eta2\n";

    ofstream prof("profiles1.csv");
    prof << "r";
    for (double T : save_points) prof << ",C_t" << T;
    prof << "\n";

    for (int it=0; it <= tp; it++) {
        double t=it*dt;

        double R_act=richardson(r, C, k);
        double R_ideal=k*Cs*(4.0*M_PI*R*R*R/3.0);
        double eta2=R_act/R_ideal;

        eta2_out << fixed << setprecision(6) << t << "," << eta2 << "\n";

        if (save_steps[it]) saved_profiles1.push_back(C);

        if (it == tp) break;

        vector<double> rhs(P);
        rhs[0]=eps/dt*C[0];
        for (int j=1; j<P-1; j++) rhs[j]=eps/dt*C[j];
        rhs[P-1]=Cs;

        C=thomas(a, b, c, rhs);
        C[P-1]=Cs; // enforce Dirichlet BC
    }

    for (int i=0; i<P; i++) {
        prof << r[i];
        for (auto &vec : saved_profiles1) prof << "," << vec[i];
        prof << "\n";
    }

    // files closed by destructor
    return 0;
}
