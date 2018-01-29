#include <numeric>
#include <vector>

typedef std::vector<unsigned short> vecui2;


void index_convert_one2d(const int I, const vecui2& dims, vecui2& indx)
{
    if(dims.size() > 11){
        std::cout << "dimension too high: "<< dims.size() << std::endl;
        exit(0);
    }
    else if(dims.size() == 1){
        indx.push_back(I);
    }
    else{
        indx.push_back( I / std::accumulate( dims.begin()+1, dims.end(), 1, std::multiplies<unsigned short>()) );
        int II = I % std::accumulate(dims.begin() + 1, dims.end(), 1, std::multiplies<unsigned short>());
        index_convert_one2d(II, vecui2(dims.begin()+1, dims.end()), indx);
    }
}

int index_convert_d2one(const vecui2& indx, const vecui2& dims)
{
    if(dims.size() > 11){
        std::cout << "dimension too high: "<< dims.size() << std::endl;
        exit(0);
    }
    else if(dims.size() == 1){
        return indx[0];
    }
    else{
        int i = indx[0] * std::accumulate( dims.begin()+1, dims.end(), 1, std::multiplies<unsigned short>() ) ;
        return i + index_convert_d2one(vecui2(indx.begin()+1, indx.end()), vecui2(dims.begin() + 1, dims.end()));
    }
}



