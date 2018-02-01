#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath> //abs()
#include <cassert> //assert()
#include <omp.h> //pragma parallel
#include <numeric> /* std::accumulate, std::inner_product */
#include <stdlib.h>  /* atoi */
#include "string_utils.cpp" // int2string
#include "spins_io.cpp" // vec read/write is 0,1 instead of -1,1 
/* void read_spin(const char * filename, std::vector<short>& vec); */
#include "indx_utils.cpp" 
/* void index_convert_one2d(const int I, const vecui2 dims, vecui2& indx)
   int index_convert_d2one(const vecui2& indx, const vecui2 dims) */



typedef std::vector<int> veci4;
typedef std::vector<short int> veci2;
typedef std::vector<unsigned short> vecui2;
typedef std::vector<int>::const_iterator veci4_iter;
typedef std::vector<short int>::const_iterator veci2_iter;
typedef std::vector<unsigned short>::const_iterator vecui2_iter;

inline int mod(int a, int b)
{
    int r = a % b;
    return r < 0 ? r + b : r;
}


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
        int e,m;
        IsingModel(std::string filename, vecui2& dims, float T, bool is_fBC = false);
        void Monte_Carlo();
        void get_energy();
        void get_m();
        void draw2D(int i);
        

};

void IsingModel::get_energy()
{
    veci4 nns;
    this->e = 0;
    for (veci4_iter iter = this->valid_pos.begin(); iter != this->valid_pos.end(); ++iter){
        this->find_nns(*iter, nns);
        for (int j = 0; j < this->dim ; ++j)
            this->e -= this->lattice[*iter] * this->lattice[nns[j]]; /* look at positive direction only */
    }
    assert(abs(this->e) <= this->valid_pos.size() * this->dim);
}

void IsingModel::get_m()
{
    this->m = 0;
    for (auto iter = this->valid_pos.begin(); iter != this->valid_pos.end(); ++iter){
        this->m += this->lattice[*iter]; /* look at positive direction only */
    }
    assert(abs(this->m) <= this->valid_pos.size());
}

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

    if(is_fBC)
        for(vecui2_iter itr = dims.begin(); itr != dims.end(); ++itr)
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
    printf(" step = %d, e = %d ", t,  this->e);
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
    // if(draw_s) ising.draw2D(0);
    std::ofstream e_file("e_dist_T");
    std::ofstream m_file("m_dist_T");
    for(int t = 0; t < tmax; ++t){
        ising.Monte_Carlo();
        ising.get_energy();
        ising.get_m();
        if(draw_s) ising.draw2D(t+1);
        e_file << ising.e << std::endl;
        m_file << ising.m << std::endl;
    }

    if(draw_s) {std::cout << std::endl; std::cout << "\033[" << dims[0]+ 4 << "B";}

    e_file.close();
    m_file.close();
}

int main(int argc, char *argv[])
{
    // set default
    vecui2 dims; /* {N,M} */
    dims.push_back(1000);dims.push_back(1000);
    float T = 2.27; //Tc = 2.269
    bool draw_s = false;
    bool is_fBC = false;
    int tmax = 2;
    std::string filename("s_pBC_L1000_1mMC");
    
  
    // argv[0]:program name, argv[1]:T, argv[2]:BC, argv[3:]:dims
    srand(222);    
    dist_gen(T, dims, tmax, filename, draw_s, is_fBC);
}


