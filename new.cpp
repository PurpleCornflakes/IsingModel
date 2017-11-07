/* 2D Ising Model */

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <cassert>
#include <numeric>      // std::accumulate, std::inner_product

typedef std::vector<int> vec;
typedef std::vector<float> vecf;

inline int mod(int a, int b);
class IsingModel
{
private:
    int N, M;
    int num_s;
    float T;
    vec s;
    float e; //current system energy
    float e_sqr;
    void find_nns(int n, int (&nns)[4]);
    void energy();
public:
    IsingModel(int N, int M);
    ~IsingModel(void);
    void Monte_Carlo(bool eval_energy);
    inline void set_T(float new_T);
    inline void get_energy(float ee[2]);
    inline void get_m(float* m);
    void draw(int i, bool eval_energy);
};

void dist_gen(float T, bool e_dist, bool m_dist, int N, int M, int tmax, int tdis)
{
    
    // if(e_dist)
        float ee[2]; // array to store energy and squared energy
        std::ofstream e_file;
        e_file.open("e_dist", std::ios::out);
    // if(m_dist)
        float m;
        std::ofstream m_file;
        m_file.open("m_dist", std::ios::out);

    IsingModel ising(N, M);
    ising.draw(0, false);
    ising.set_T(T);
    for(int t = 0; t < tdis; ++t){
        ising.Monte_Carlo(false);
        ising.draw(t+1, false);}

    for(int t = 0; t < tmax - tdis; ++t){
        ising.Monte_Carlo(true);
        ising.draw(t+tdis+1, false);
        ising.get_energy(ee);
        ising.get_m(&m);
        if(e_dist) 
            e_file << ee[0] << std::endl;
        if(m_dist)  
            m_file << m << std::endl;}
        e_file.close();
        m_file.close();
}

void e_T_m_T_gen(bool e_T, bool m_T, int N, int M,int tmax, int tdis){
    int ave_n = tmax - tdis;
    vecf T;
    for (float tt = 2.0; tt < 3.5 ; tt += 0.2)
        T.push_back(tt);
    // if(e_T){
        float c;
        float ee[2]; 
        vecf e_sum, e2_sum;
        float e_ave, e2_ave; //average energy, average squared energy.
        std::ofstream e_T_file;
        e_T_file.open("e_T", std::ios::out);
    // if(m_T){
        float chi, m;
        vecf m_sum, m2_sum;
        float m_ave, m2_ave; //average tol_spin, average squared tol_spin.
        std::ofstream m_T_file;
        m_T_file.open("m_T", std::ios::out);

    IsingModel ising(N, M);
    ising.draw(0, false);
    for(int i = 0; i < T.size(); i++){
        ising.set_T(T[i]);
        if(e_T){
            e_sum.push_back(0.);
            e2_sum.push_back(0.);}
        if(m_T){
            m_sum.push_back(0);
            m2_sum.push_back(0);}
        for(int t = 0; t < tdis; ++t){
            ising.Monte_Carlo(false);
            ising.draw(t+1, false);}
        for(int t = 0; t < tmax - tdis; ++t){
            ising.Monte_Carlo(true);
            ising.draw(t+tdis+1, false);
            ising.get_energy(ee);
            ising.get_m(&m);
            if(e_T) {
                e_sum[i] += ee[0];
                e2_sum[i] += ee[1];}
            if(m_T){
                m_sum[i] += m;
                m2_sum[i] += m*m;}
            }
        if(e_T){
            e_ave = e_sum[i]/ave_n; e2_ave = e2_sum[i]/ave_n;
            c = (e2_ave - e_ave * e_ave)/T[i]/T[i];
            e_T_file << e_ave << "," << c << "\n";}
        if(m_T){
            m_ave = m_sum[i]/ave_n; m2_ave = m2_sum[i]/ave_n;
            chi = m2_ave - m_ave*m_ave;
            m_T_file << m_ave << "," << chi << "\n";}
        }
        e_T_file.close();
        m_T_file.close();
}

int main()
{
    int N = 32, M = 32;
    int tmax = 5000;
    int tdis = 1000;
    bool e_dist = true;
    bool m_dist = false; 
    bool e_T = false;
    bool m_T = true;
    float T = 1.5;

    // dist_gen(T, e_dist, m_dist, N, M, tmax, tdis);
    e_T_m_T_gen(e_T, m_T, N, M, tmax, tdis);

    std::cout << "\033[" << N+1 << "B" <<std::endl;

}




// constructor
IsingModel::IsingModel(int N, int M)
{
    this->N = N;
    this->M = M;
    this->num_s = N*M;
    srand(222);
    for(int i = 0; i < N*M; ++i)
        this->s.push_back( (int)(2*(std::rand()%2 - 0.5)) );
}

IsingModel::~IsingModel()
{
    ;
}


void IsingModel::Monte_Carlo(bool eval_energy = false)
{
    int nns[4];
    int e = 0;
    int n;
    for(int k = 0; k < this->num_s; ++k){
        n = std::rand()%this->num_s;
        this->find_nns(n, nns);
        e = 0;
        for(int j = 0; j < 4; ++j)
            e += this->s[nns[j]];
        e *= 2 * this->s[n];
        if(e <= 0 || std::rand()/(float)RAND_MAX < std::exp(-e/this->T))
            this->s[n] *= -1;
    }
    // evaluate energy
    if (eval_energy) this->energy();
}

void IsingModel::find_nns(int n, int (&nns)[4])
{

    int i, j, u_i, d_i, l_j, r_j;

    i = n/this->M; j = n%this->M;
    u_i = mod(i-1, this->N); 
    d_i = mod(i+1, this->N);
    l_j = mod(j-1, this->M); 
    r_j = mod(j+1, this->M);
    
    nns[0] = u_i*this->N + j; 
    nns[1] = d_i*this->N + j;
    nns[2] = i*this->M + l_j; 
    nns[3] = i*this->M + r_j;
}

void IsingModel::energy()
{
    int nns[4];
    this->e = 0;
    for (int i = 0; i < this->num_s; ++i){
        this->find_nns(i, nns);
        this->e -= this->s[i] * this->s[nns[1]];
        this->e -= this->s[i] * this->s[nns[3]];
    }
    assert(abs(this->e) <= this->num_s * 2 + 1e-5);
    this->e_sqr = this->e * this->e;
}

inline void IsingModel::get_energy(float ee[2])
{
    ee[0] = this->e;
    ee[1] = this->e_sqr;
}

inline void IsingModel::get_m(float *m)
{
    *m = std::accumulate(this->s.begin(), this->s.end(), 0);
}

inline void IsingModel::set_T(float new_T)
{
    this->T = new_T;
}

void IsingModel::draw(int i, bool eval_energy = false)
{
    std::string str = "\n";
    for (int i = 0; i < this->num_s ; ++i){
        str += this->s[i] == 1 ? "\033[01;31m O\033[00m" : "\033[01;32m O\033[00m";
        if ((i+1) % this->M == 0) str += "\n";
    }
    std::cout << str;
    printf(" step = %d, T = %.3f ", i,  this->T);
    if (eval_energy)
        printf(", E = %.0f", this->e );
    std::fflush(stdout);
    std::cout << "\033[" << this->N+1 << "A";
}

inline int mod(int a, int b)
{
    int r = a % b;
    return r < 0 ? r + b : r;
}


