#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

void write_spin(const char * filename, const std::vector<short>& vec);
void read_spin(const char * filename, std::vector<short>& vec);


// int main(){
//     std::vector<short> vec;
//     std::vector<short> vec2;
//     for (int i = 0; i < 1000000; ++i){
//         vec.push_back(i%2 - (i+1)%2);
//     }

//     write_spin("test.bin", vec);
//     read_spin("test.bin", vec2);
//     for (int i = 0; i < 100; ++i){
//         std::cout << vec2[i] << ",";
//     }
// }

void write_spin(const char * filename, const std::vector<short>& vec)
{
    std::ofstream outfile(filename);
    for (auto itr = vec.begin(); itr!=vec.end(); ++itr){
        outfile << *itr;
    }
    outfile.close();
}

void read_spin(const char * filename, std::vector<short>& vec)
{
    std::ifstream infile(filename);
    if(infile){
        infile.seekg(0, infile.end);
        int size = infile.tellg();
        infile.seekg(0, infile.beg);
        char * buffer = new char [size];
        // std::cout << "file size = " << size << "bytes, reading ..." << std::endl;
        infile.read(buffer, size);
        infile.close();

        vec.clear();
        for (int i = 0; i< size; ++i){
            vec.push_back(buffer[i] - '0');
        }
    }
    else{
        std::cout << "cannot open file" << std::endl;
    }
}