#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <cassert>
#include <omp.h>
#include <numeric>
/* std::accumulate, std::inner_product */
#include <stdlib.h>  /* atoi */
#include "string_utils.cpp"
#include "spins_io.cpp" // vec read/write is 0,1 instead of -1,1 
/* void read_spin(const char * filename, std::vector<short>& vec); */
#include "indx_utils.cpp" 
/* void index_convert_one2d(const int I, const vecui2 dims, vecui2& indx)
   int index_convert_d2one(const vecui2& indx, const vecui2 dims) */



typedef std::vector<int> veci4;
typedef std::vector<short int> veci2;
typedef std::vector<unsigned short> vecui2;


inline int mod(int a, int b)
{
    int r = a % b;
    return r < 0 ? r + b : r;
}


// void IsingModel::energy()
// {
//     int nns[4];
//     this->e = 0;
//     for (vecl_iter iter = this->pos_valid.begin(); iter != this->pos_valid.end(); ++iter){
//         this->find_nns(*iter, nns);
//         this->e -= this->s[*iter] * this->s[nns[1]]; /* down */
//         this->e -= this->s[*iter] * this->s[nns[3]]; /* right */
//     }
//     assert(abs(this->e) <= this->num_s * 2 + 1e-5);
//     this->e_sqr = this->e * this->e;
// }

class IsingModel
{
    private:
        inline void find_nns(int I, veci4& nns);
        float T;
        int num_lattice;
        int dim;
        vecui2 dims; // dimensions
        veci4 valid_pos;
        veci2 lattice;

    public: 
        IsingModel(std::string filename, vecui2& dims, float T, bool is_fBC = false);
        void Monte_Carlo();
        void draw2D(int i);
};



IsingModel::IsingModel(std::string filename, vecui2& dims, float T, bool is_fBC)
{
    this->T = T;
    // read s_file(0,1) & convert to lattice
    veci2 s_vec10;
    const char* c_filename = filename.c_str();
    read_spin(c_filename,s_vec10);

    this->num_lattice = s_vec10.size();
    this->lattice.resize(this->num_lattice);

    for(int i = 0; i < this->num_lattice; ++i){
        this->lattice[i] = (s_vec10[i] == 1 ? 1 : -1);}

    // for(auto iter = this->lattice.begin(); iter != this->lattice.end(); ++iter)
    //     std::cout << *iter << std::endl;

    if(is_fBC)
        for(auto itr = dims.begin(); itr != dims.end(); ++itr)
            this->dims.push_back(*itr + 2);
    else 
        this->dims = dims;
    this->dim = this->dims.size();

    this->valid_pos.clear();
    if(is_fBC)   
        for(int I = 0; I < this->num_lattice; ++I){
            vecui2 indx; indx.clear(); 
            index_convert_one2d(I, this->dims, indx);
            //I is valid if not on boundary
            bool is_valid = true;
            for(int d = 0; d < this->dim; ++d)
                if(indx[d] == 0 || indx[d] == this->dims[d]-1){
                    is_valid = false;
                    this->lattice[I] = 0;
                    break;
                }
            if(is_valid)
                this->valid_pos.push_back(I);
        }
    else
        for(int I = 0; I < this->num_lattice; ++I)
            this->valid_pos.push_back(I);
}


inline void IsingModel::find_nns(int I, veci4& nns)
{
    nns.clear();
    vecui2 indx;
    index_convert_one2d(I, this->dims, indx);

    vecui2 nn_indx(indx);
    for(int d = 0; d < this->dim; ++d){
        nn_indx[d] = mod(indx[d]+1, this->dims[d]);
        nns.push_back(index_convert_d2one(nn_indx, this->dims));
        nn_indx[d] = indx[d];
    }
    for(int d = 0; d < this->dim; ++d){
        nn_indx[d] = mod(indx[d]-1, this->dims[d]);
        nns.push_back(index_convert_d2one(nn_indx, this->dims));
        nn_indx[d] = indx[d];
    }

}

void IsingModel::draw2D(int t)
{
    int N = this->dims[0]; int M= this->dims[1];
    std::string str = "\n";
    for (int i = 0; i < this->num_lattice ; ++i){
        str += this->lattice[i] == 1 ? "\033[01;31m O\033[00m" : "\033[01;32m O\033[00m";
        if ((i+1) % M == 0) str += "\n";
    }
    std::cout << str;
    printf(" step = %d, T = %.3f ", t,  this->T);
    std::fflush(stdout);
    std::cout << "\033[" << N+1 << "A";
}

void IsingModel::Monte_Carlo()
{
    veci4 nns;
    int e, I;
    for(int k = 0; k < this->valid_pos.size(); ++k){ 
        I = this->valid_pos[std::rand()%this->valid_pos.size()];
        this->find_nns(I, nns);

        for(auto iter= nns.begin(); iter != nns.end(); ++iter)
            std::cout << I << ":" << *iter << std::endl;

        e = 0;
        for(int j = 0; j < nns.size(); ++j)
            e += this->lattice[nns[j]];
        e *= 2 * this->lattice[I];

        if(e <= 0 || std::rand()/(float)RAND_MAX < std::exp(-e/this->T))
            this->lattice[I] *= -1;
    }
}



void dist_gen(float T, vecui2& dims, int tmax, std::string filename, bool draw_s, bool is_fBC)
{
    IsingModel ising(filename, dims, T, is_fBC);
    if(draw_s) ising.draw2D(0);
    if(draw_s) std::cout << "\033[" << dims[0]+3 << "B" <<std::endl;
    // float e; // array to store energy and squared energy
    // std::ofstream e_file;
    // e_file.open("e_dist_T15", std::ios::out);
    // float m;
    // std::ofstream m_file;
    // m_file.open("m_dist_T15", std::ios::out);
    // for(int t = 0; t < tmax; ++t){
    //     ising.Monte_Carlo();
    //     ising.draw2D(t+1);
    //     ising.get_energy(e);
    //     ising.get_m(&m);
    //     e_file << e << std::endl;
    //     m_file << m << std::endl;
    // }
    // e_file.close();
    // m_file.close();
}

int main(int argc, char *argv[])
{
    // set default
    vecui2 dims({5,6}); /* {N,M} */
    float T = 2.27; //Tc = 2.269
    bool draw_s = true;
    bool is_fBC = false;
    int tmax = 100;
    std::string filename("ss_T2.25_pBC_d5_d6_30MC");
    
  
    // argv[0]:program name, argv[1]:T, argv[2]:BC, argv[3:]:dims
    srand(222);    
    dist_gen(T, dims, tmax, filename, draw_s, is_fBC);
}


