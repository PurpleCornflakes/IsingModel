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
/* void write_spin(const char * filename, const std::vector<short>& vec) */
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
        IsingModel(vecui2& dims, float T, bool is_fBC = false);
        void Monte_Carlo();
        void draw2D(int i);
        void Write_lattice(std::string filename);
};



IsingModel::IsingModel(vecui2& dims, float T, bool is_fBC)
{
    if(is_fBC)
        for(auto itr = dims.begin(); itr != dims.end(); ++itr)
            this->dims.push_back(*itr + 2);
    else 
        this->dims = dims;

    this->dim = this->dims.size();
    this->num_lattice =std::accumulate(this->dims.begin(), this->dims.end(), 1, std::multiplies<unsigned short>());
    this->valid_pos.clear();
    if(is_fBC){
        #pragma omp parallel
        {
            #pragma omp for ordered
            for(int I = 0; I < this->num_lattice; ++I){
                vecui2 indx; indx.clear(); 
                index_convert_one2d(I, this->dims, indx);
                //I is valid if not on boundary
                bool is_valid = true;
                for(int d = 0; d < this->dim; ++d)
                    if(indx[d] == 0 || indx[d] == this->dims[d]-1){
                        is_valid = false; 
                        break;
                    }
                if(is_valid){
                    #pragma omp ordered
                    this->valid_pos.push_back(I);
                }
            }
        }
    }
    else
        for(int I = 0; I < this->num_lattice; ++I)
            this->valid_pos.push_back(I);

    this->T = T;
    this->lattice.resize(this->num_lattice, 0); // init all to 0
    for(auto iter = this->valid_pos.begin(); iter != this->valid_pos.end(); ++iter)
        this->lattice[*iter] = 1;
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

        e = 0;
        for(int j = 0; j < nns.size(); ++j)
            e += this->lattice[nns[j]];
        e *= 2 * this->lattice[I];

        if(e <= 0 || std::rand()/(float)RAND_MAX < std::exp(-e/this->T))
            this->lattice[I] *= -1;
    }
}

void IsingModel::Write_lattice(std::string filename)
{
    // write lattice(-1,1) into file(0,1)
    veci2 s_vec10;
    for(auto iter = this->lattice.begin(); iter != this->lattice.end(); ++iter)
        s_vec10.push_back(*iter == 1 ? 1 : 0);
    const char* c_filename = filename.c_str();
    write_spin(c_filename,s_vec10);
}

void thermalize(float T, vecui2& dims, int tdis, std::string filename, bool draw_s, bool is_fBC)
{
    
    IsingModel ising(dims, T, is_fBC);

    for(int t = 0; t < tdis; ++t){
        ising.Monte_Carlo();
        if(draw_s) ising.draw2D(t+1);
    }
    if(draw_s) std::cout << "\033[" << dims[0]+3 << "B" <<std::endl;
    ising.Write_lattice(filename);
}


int main(int argc, char *argv[])
{
    // set default
    vecui2 dims({2,3}); /* {N,M} */
    float T = 2.27; //Tc = 2.269
    bool draw_s = false;
    bool is_fBC = true;
    int tdis = 3;
    std::string filename("ss");
    
  
    // argv[0]:program name, argv[1]:T, argv[2]:BC, argv[3:]:dims
    if (argc > 1){
        T = atof(argv[1]);
        std::string str_T(argv[1]), str_BC(argv[2]), str_L(argv[3]);
        if(str_BC == "fBC") 
            is_fBC = true;
        else
            is_fBC = false;

        int dim = argc-3;
        dims.resize(dim);
        filename += "_T" + str_T + "_" + str_BC + "_";
        for(int i = 0; i < dim; ++i){
            dims[i] = atoi(argv[i+3]);
            tdis *= dims[i];
            filename += "d" + mine::int2string(dims[i]) + "_";
        }
        filename += mine::int2string(tdis)+ "MC";
    }
    srand(222);    
    thermalize(T, dims, tdis, filename, draw_s, is_fBC);
}
